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

Particle ParticleFilter::AddGaussianNoise(Particle& particle, const double std_pos[]) {
    normal_distribution<double> dist_x(particle.x, std_pos[0]);
    normal_distribution<double> dist_y(particle.y, std_pos[1]);
    normal_distribution<double> dist_theta(particle.theta, std_pos[2]);

    particle.x = dist_x(this->gen);
    particle.y = dist_y(this->gen);
    particle.theta = dist_theta(this->gen);
    return particle;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    if(!this->is_initialized) {
        this->num_particles = 10;

        for(int i = 0; i < this->num_particles; i++) {
            this->weights.push_back(1./this->num_particles);
            Particle particle {
                i, // id,
                x, //x
                y, //y
                theta, //theta
                1./this->num_particles // weight
            };

            this->particles.push_back(ParticleFilter::AddGaussianNoise(particle, std));
        }

        this->is_initialized = true;
    }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    for(Particle& particle : this->particles) {
        if(abs(yaw_rate) < 1e-6) {
            particle.x += delta_t * velocity * cos(particle.theta);
            particle.y += delta_t * velocity * sin(particle.theta);
        } else {
            particle.x += velocity/yaw_rate * (sin(particle.theta + delta_t * yaw_rate) - sin(particle.theta));
            particle.y += velocity/yaw_rate * (cos(particle.theta) - cos(particle.theta + delta_t * yaw_rate));
        }
        particle.theta += yaw_rate * delta_t;
        particle = ParticleFilter::AddGaussianNoise(particle, std_pos);
    }
}

std::vector<LandmarkObs> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    std::vector<LandmarkObs> result;
    for(LandmarkObs& observation : observations) {
        double minDist = numeric_limits<float>().max();
        const LandmarkObs* closest = nullptr;
        for(const LandmarkObs& pred : predicted) {
            double distance = dist(observation.x, observation.y, pred.x, pred.y);
            if(distance < minDist) {
                minDist = distance;
                closest = &pred;
            }
        }

        if(closest != nullptr) {
            result.push_back(*closest);
        }
    }

    return result;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
        const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

    double gauss_norm= 1./(2 * M_PI * std_landmark[0] * std_landmark[1]);

    for(int i = 0;i < num_particles; i++) {
        Particle& particle = this->particles[i];

        vector<LandmarkObs> predicted;
        for(const auto& landmark : map_landmarks.landmark_list) {
            if(fabs(landmark.x_f - particle.x) < sensor_range && fabs(landmark.y_f - particle.y) < sensor_range) {
                LandmarkObs landmarkObs {
                    landmark.id_i, //id
                    landmark.x_f,  //x
                    landmark.y_f   //y
                };
                predicted.push_back(landmarkObs);
            }
        }

        vector<LandmarkObs> observationsInMap;
        for(const LandmarkObs& observation : observations) {
            LandmarkObs observationInMap {
                observation.id, //id
                particle.x + cos(particle.theta) * observation.x - sin(particle.theta) * observation.y, //x
                particle.y + sin(particle.theta) * observation.x + cos(particle.theta) * observation.y  //y
            };
            observationsInMap.push_back(observationInMap);
        }

        std::vector<LandmarkObs> associatedObservation = dataAssociation(predicted, observationsInMap);
        vector<int> associations;
        vector<double> sense_x, sense_y;
        double weight = 1.0;
        for(size_t j = 0; j < observationsInMap.size(); j++) {
            LandmarkObs landmarkObs = associatedObservation[j];
            LandmarkObs observation = observationsInMap[j];

            associations.push_back(landmarkObs.id);
            sense_x.push_back(landmarkObs.x);
            sense_y.push_back(landmarkObs.y);
            
            double exponent = pow(landmarkObs.x - observation.x, 2)/(2 * pow(std_landmark[0], 2)) + pow(landmarkObs.y - observation.y, 2)/(2 * pow(std_landmark[1], 2));
            weight *= gauss_norm * exp(-exponent);
        }

        this->SetAssociations(particle, associations, sense_x, sense_y);
        this->weights[i] = weight;
    }
}

void ParticleFilter::resample() {
    vector<Particle> result;
    discrete_distribution<> distribution(begin(this->weights), end(this->weights));

    for(int i = 0; i < this->num_particles; i++) {
        result.push_back(this->particles[distribution(this->gen)]);
    }

    this->particles = result;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
