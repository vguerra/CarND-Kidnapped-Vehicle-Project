/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 100;

  particles = std::vector<Particle>(num_particles);

  std::default_random_engine generator;

  std::normal_distribution<double> x_distribution(x, std[0]);
  std::normal_distribution<double> y_distribution(y, std[1]);
  std::normal_distribution<double> theta_distribution(theta, std[2]);

  int id = 0;
  for (auto& particle : particles)
  {
    particle.id = id;
    particle.x = x_distribution(generator);
    particle.y = y_distribution(generator);
    particle.theta = theta_distribution(generator);
    particle.weight = 1.0;
    id += 1;
  }

  weights = std::vector<double>(num_particles, 1.0);

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  std::default_random_engine generator;
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  if (yaw_rate == 0) {
    // bicycle model with zero yaw_rate
    double factor = velocity * delta_t;

    for (auto& particle : particles) {
      particle.x += factor * cos(particle.theta);
      particle.y += factor * sin(particle.theta);

      particle.x = std::normal_distribution<double>(particle.x, std_x)(generator);
      particle.y = std::normal_distribution<double>(particle.y, std_y)(generator);
      particle.theta = std::normal_distribution<double>(particle.theta, std_theta)(generator);
    }

  } else {
    // bicycle model
    double factor = velocity/yaw_rate;
    double delta_theta = yaw_rate*delta_t;

    for (auto& particle : particles) {
      double new_theta = particle.theta + delta_theta;
      particle.x += factor*(sin(new_theta) - sin(particle.theta));
      particle.y += factor*(cos(particle.theta) - cos(new_theta));
      particle.theta = new_theta;

      particle.x = std::normal_distribution<double>(particle.x, std_x)(generator);
      particle.y = std::normal_distribution<double>(particle.y, std_y)(generator);
      particle.theta = std::normal_distribution<double>(particle.theta, std_theta)(generator);
    }
  }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  // Constants through-out the process of computing weights
  double std_x = std_landmark[0];
  double std_x_2 = std_x*std_x;
  double std_y = std_landmark[1];
  double std_y_2 = std_y*std_y;
  double multivariate_factor = 1/(2*M_PI*std_x*std_y);


  for (auto& particle : particles) {

//    // Landmarks in sensor range
//    std::vector<Map::single_landmark_s> landmarks_within_range;
//    landmarks_within_range.reserve(map_landmarks.landmark_list.size());
//
//    std::copy_if(map_landmarks.landmark_list.begin(),
//                 map_landmarks.landmark_list.end(),
//                 std::back_inserter(landmarks_within_range),
//                 [&particle](const Map::single_landmark_s& landmark) {
//                   return
//
//                 });

    // Transform observations first to map coordinates
    std::vector<LandmarkObs> transformed_observations(observations.size());

    std::transform(observations.begin(),
                   observations.end(),
                   transformed_observations.begin(),
                   [&particle](LandmarkObs observation) {

                     LandmarkObs transformed_observation;

                     transformed_observation.x = particle.x +
                     observation.x*cos(particle.theta) -
                     observation.y*sin(particle.theta);

                     transformed_observation.y = particle.y +
                     observation.x*sin(particle.theta) +
                     observation.y*cos(particle.theta);

                     return transformed_observation;
                   });

    vector<int> associations(transformed_observations.size());
    vector<double> sense_x(transformed_observations.size());
    vector<double> sense_y(transformed_observations.size());

    // Look for closest Map landmark for each transformed observation
    for (size_t i = 0; i < transformed_observations.size(); ++i) {

      double min_distance = std::numeric_limits<double>::max();;
      int best_landmark_id = 0;
      double best_x = 0.0;
      double best_y = 0.0;

      for (auto landmark : map_landmarks.landmark_list) {
        double distance = dist(transformed_observations[i].x, transformed_observations[i].y, landmark.x_f, landmark.y_f);
        if (distance < min_distance) {
          min_distance = distance;
          best_landmark_id = landmark.id_i;
          best_x = landmark.x_f;
          best_y = landmark.y_f;
        }
      }

      associations[i] = best_landmark_id;
      sense_x[i] = best_x;
      sense_y[i] = best_y;
    }

    particle = SetAssociations(particle, associations, sense_x, sense_y);

    double particle_weight = 1.0;

    for (size_t i = 0; i < particle.associations.size(); ++i) {
      LandmarkObs observation = transformed_observations[i];
      double exponent = pow(observation.x - sense_x[i], 2.0)/(2*std_x_2) +
          pow(observation.y - sense_y[i], 2.0)/(2*std_y_2);
      double res = exp(-exponent);
      particle_weight *= multivariate_factor*res;
    }

    particle.weight = particle_weight;
    weights[particle.id] = particle_weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::default_random_engine generator;
  std::discrete_distribution<double> dist(weights.begin(), weights.end());

  std::vector<Particle> sampled_particles(particles.size());

  for (size_t p = 0; p < particles.size(); ++p) {
    int drawn_index = dist(generator);
    sampled_particles[p] = particles[drawn_index];
    sampled_particles[p].id = p;
  }

  particles = sampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
