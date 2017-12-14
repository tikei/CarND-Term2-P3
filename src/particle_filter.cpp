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
    num_particles = 70;
    weights.resize(num_particles);

    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];

    // random engine; random device rd is initialized in header file
    std::mt19937 rand_engine(rd());
    std::normal_distribution<double> dist_x(x, std_x);
    std::normal_distribution<double> dist_y(y, std_y);
    std::normal_distribution<double> dist_theta(theta, std_theta);

    // Initialize all particles, drawing from a Normal distribution
    for (int i = 0; i < num_particles; ++i) {
        double sample_x, sample_y, sample_theta;

        sample_x = dist_x(rand_engine);
        sample_y = dist_y(rand_engine);
        sample_theta = dist_theta(rand_engine);

        Particle particle;

        particle.id = i;
        particle.x = sample_x;
        particle.y = sample_y;
        particle.theta = sample_theta;
        particle.weight = 1.0f;

        particles.push_back(particle);
        weights[i] = particle.weight;

    }// end for

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // double std_a = std_pos[0]; // std deviation of velocity(control noise)
    // double std_yawd = std_pos[1]; // standard deviation of yaw rate;

    // std::normal_distribution<double> dist_a(velocity, std_pos[0]);
    // std::normal_distribution<double> dist_yaw(yaw_rate, std_pos[1]);

    // add random noise
    std::mt19937 rand_engine(rd());
    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0, std_pos[2]);

    double scale = (velocity / yaw_rate);
    for (int i = 0; i < particles.size(); i++){
        double theta = particles[i].theta;

        // avoid division by zero
        if (std::fabs(yaw_rate) < 0.0001){
            particles[i].x += velocity * delta_t * std::cos(theta);
            particles[i].y += velocity * delta_t * std::sin(theta);
        } else {
            particles[i].x += scale * (std::sin( theta + yaw_rate * delta_t) - std::sin(theta));
            particles[i].y += scale * (std::cos(theta) - std::cos( theta + yaw_rate * delta_t));
            particles[i].theta += yaw_rate * delta_t;

        }//if else

        // update with predicted values, adding noise
        particles[i].x += dist_x(rand_engine);
        particles[i].y += dist_y(rand_engine);
        particles[i].theta += dist_theta(rand_engine);
    }//for
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

    // STEP 1:
    // find each particle find the landmarks within sensor_range
    // associated with the particle
    for (int i = 0; i < particles.size(); i++){

        // particle coordinates
        double x_part = particles[i].x;
        double y_part = particles[i].y;
        double theta_part = particles[i].theta;

        double total_prob = 1.0;

        for (int j = 0; j < observations.size(); j++){
            // get the observation coordinates
            double x_obs = observations[j].x;
            double y_obs = observations[j].y;

            // STEP 2:
            // rotate and translate the sensor observation coordinates in map coordinates
            double x_map = x_part + (std::cos(theta_part) * x_obs) - (std::sin(theta_part) * y_obs);
            double y_map = y_part + (std::sin(theta_part) * x_obs) + (std::cos(theta_part) * y_obs);

            // STEP 3
            // find nearest landmark to the sensor observation
            double min_distance = std::numeric_limits<double>::max();
            int landmark_idx;

            for (int k = 0; k < map_landmarks.landmark_list.size(); k++){
                // filter: consider only landmarks in range:
                bool in_x_range = false;
                bool in_y_range = false;
                Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
                if (x_part >= (landmark.x_f - sensor_range) &&
                        x_part <= (landmark.x_f + sensor_range)){
                    in_x_range = true;
                }//if
                if (y_part >= (landmark.y_f - sensor_range) &&
                        y_part <= (landmark.y_f + sensor_range)){
                    in_y_range = true;
                }//if
                if (in_x_range && in_y_range){
                    // calculate minimum distance observation <> landmark

                    // possibility to calculate L-1 norm as
                    // L-2 (Euclidean) more expensive due to sqrt
                    // double distance = std::abs(x_map - landmark.x_f) +
                    //     std::abs(y_map - landmark.y_f);
                    double distance = dist(x_map, y_map, landmark.x_f, landmark.y_f);

                    if (distance < min_distance){
                        min_distance = distance;
                        // the k index of observations represents the
                        // closesest landmark to the sensor observation
                        landmark_idx = k;
                    }//if
                }//if in range
            }//end for k landmarks

            // mu_x , mu_y = nearest landmark x and y
            double mu_x = map_landmarks.landmark_list[landmark_idx].x_f;
            double mu_y = map_landmarks.landmark_list[landmark_idx].y_f;

            // calculate probability for this association observation<>landmark
            double prob = MultiVariateNormalPdf(x_map, y_map, mu_x, mu_y, std_landmark[0], std_landmark[1]);

            // total probability of the particle based on the probabilities
            // of each association observation <> landmark:
            total_prob *= prob;

        } //end for j observations
        // update particle's weight
        particles[i].weight = total_prob;

        // update all Particle filter weights
        weights[i] = total_prob;
    }// end for i
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // temp particles vector
    std::vector<Particle> resampled_particles;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> dist(weights.begin(), weights.end());
    for (int i = 0; i < num_particles; i++){
        resampled_particles.push_back(particles[dist(gen)]);
    }//end for

    particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
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
