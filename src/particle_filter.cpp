/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

//static default_random_engine gen;

using std::string;
using std::vector;

//using namespace std;

using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  is_initialized=true; 
  
  num_particles = 10;  // TODO: Set the number of particles
  
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  // TODO: Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2]; 
  
  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x, y and theta
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);
  
  //Init Particles
  for (int i=0; i<num_particles; i++)
  {
    Particle particle;
    //particle.id=i;
    //particle.x=x;
    //particle.y=y;
    //particle.theta=theta;
    particle.weight=1.0;
    
    // Noise
    // where "gen" is the random engine
    particle.x=dist_x(gen);
    particle.y=dist_y(gen);
    particle.theta=dist_theta(gen);
    
    particles.push_back(particle);
  }
 //is_initialized = true;
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // Extracting standard deviations
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0, std_x);
  std::normal_distribution<double> dist_y(0, std_y);
  std::normal_distribution<double> dist_theta(0, std_theta);
  
  for (int i=0; i<num_particles; i++)
  {
    // calculate new state
    if (fabs(yaw_rate) < 0.00001) // when no change in theta
    {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } 
    else 
    {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)); // Lesson 5, Section 8
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)); // Lesson 5, Section 8
      particles[i].theta += yaw_rate * delta_t;  // Lesson 5, Section 8
    }
      //Noise
      particles[i].x += dist_x(gen);
      particles[i].y += dist_y(gen);
      particles[i].theta += dist_theta(gen);
    
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) 
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
    for (unsigned int i=0; i<observations.size(); i++) 
    {
      double min_dist=std::numeric_limits<double>::max();
      int map_id=-1;
      
      for(unsigned int j=0; j<predicted.size(); j++) 
      {
        double xDistance=observations[i].x-predicted[j].x;
        double yDistance=observations[i].y-predicted[j].y;
      
      	//double distance = (observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
        double distance = (xDistance*xDistance) + (yDistance*yDistance);
      
      	if (distance < min_dist)
      	{
          min_dist=distance;
          map_id=predicted[j].id;
        }
      }
      observations[i].id = map_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double predictions_x;
  double predictions_y;
  
  for (int i=0; i<num_particles; i++)
  {
    double particle_x=particles[i].x;
    double particle_y=particles[i].y;
    double particle_theta=particles[i].theta;
    
     // Vector to hold the map landmark locations predicted to be within sensor range of the particle
    vector<LandmarkObs> predictions;
    
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++)
    {
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;
      
      // Considering landmarks within the sensor range in rectangular region (not in circle)
      if (fabs(lm_x - particle_x) <= sensor_range && fabs(lm_y - particle_y) <= sensor_range)
      {
        predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
      }
    }
    
    //Converting (rotation and translation) landmarks observations
    vector<LandmarkObs> transformed_os;
    for (unsigned int j=0; j<observations.size(); j++)
    {
      double t_x=cos(particle_theta)*observations[j].x-sin(particle_theta)*observations[j].y+particle_x;
      double t_y=cos(particle_theta)*observations[j].x+sin(particle_theta)*observations[j].y+particle_y;
      transformed_os.push_back(LandmarkObs{observations[j].id, t_x, t_y});
    }
    dataAssociation(predictions, transformed_os);
    
    //Calculating the Particle's Final Weight; Lesson 5, Section 19
    particles[i].weight=1.0;
    weights[i]=1U;
    for (unsigned int j=0; j<transformed_os.size(); j++)
    {
      double observation_x=transformed_os[j].x;
      double observation_y=transformed_os[j].y;
      
      int associated_prediction=transformed_os[j].id;
      for (unsigned int k=0; k<predictions.size(); k++)
      {
        if (predictions[k].id==associated_prediction)
        {
          predictions_x=predictions[k].x;
          predictions_y=predictions[k].y;
          break;
        }
      }
      
      // calculate weight for this observation with multivariate Gaussian
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double dX = observation_x - predictions_x;
      double dY = observation_y - predictions_y;
      double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp(-(dX*dX/(2*s_x*s_x) + (dY*dY/(2*s_y*s_y))));
      // product of this obersvation weight with total observations weight
      particles[i].weight *=obs_w;
    }
  }      

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // get weights
  std::default_random_engine gen;
  vector<double> weights;
  double max_weight=std::numeric_limits<double>::min();
  for (int i=0; i<num_particles; i++)
  {
    weights.push_back(particles[i].weight);
    if(max_weight<particles[i].weight)
    {
      max_weight = particles[i].weight;
    }
  }
  
  // creating distributions
  std::uniform_int_distribution<int> uniintdist(0, num_particles-1);  
  std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

  int index=uniintdist(gen);
  double randomweight =unirealdist(gen);
  double beta=0.0;
  vector<Particle> resampledParticles;
  
  for (int i=0; i<num_particles; i++)
  {
    beta+=randomweight*2.0;
    while (beta> weights[index])
    {
      beta-=weights[index];
      index =(index+1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }
  particles=resampledParticles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}