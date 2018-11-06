#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

std::vector<double> JMT(std::vector<double> start, std::vector<double> end, double T)
{
  Eigen::Matrix3d A123;
  Eigen::Vector3d A345;
  Eigen::Vector3d b;
  double eq0;
  double eq1;
  double eq2;


  eq0 = end[0] - (start[0] + start[1]*T + 0.5*start[2]*T*T);
  eq1 = end[1] - (start[1] + start[2]*T);
  eq2 = end[2] - (start[2]);

  A123 <<   pow(T,3),    pow(T,4),  1*pow(T,5),
          3*pow(T,2),  4*pow(T,3),  5*pow(T,4),
          6*pow(T,1), 12*pow(T,2), 20*pow(T,3);
  
  b << eq0, eq1, eq2;
  
  // std::cout << A123 << std::endl;
  // std::cout << b << std::endl;
  
  
  A345 = A123.colPivHouseholderQr().solve(b);
  
  // std::cout << A345 << std::endl;
  
  
  return {start[0],start[1],0.5*start[2],A345(0),A345(1),A345(2)};
}

// 
bool CheckLaneBoundary(double d, int lane)
{
  bool flag = false;

  // Check that the car's d value is in the target lane.
  if(d < (2+4*lane+2) && d > (2+4*lane-2) && lane >=0 && lane <= 2)
  {
    flag = true;
  }

  return flag;
}

double SetEgoVehicleSpeed(double current_speed, double target_speed, double max_speed)
{
  double value = 0;
  // Try to match the velocity of the car in front of the ego vehicle.
  // This implementation is not perfect as it will overshoot the target by a few mph.
  if(target_speed < current_speed)
  {
    value -= 0.324;
  }
  else if(value < max_speed)
  {
    value += 0.224;
  }

  return value;
}

#define MAX_SPEED 49.5;

enum VehicleStates
{
  KeepState = 0,
  HoldLaneState,
  ChangeLanesState
};

enum SpeedStates
{
  IncreaseSpeed = 0,
  ReduceSpeed,
  MatchLaneSpeed,
  MaintainSpeed
};

// Reference velocity
double ref_vel = 0.0;

//Global loop counter
int loop_count = 0;

// Starting lane
int current_lane = 1;

// Starting states
VehicleStates state = KeepState;
SpeedStates speed_state = IncreaseSpeed;



int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

            // Sensor Fusion Parameter Order
            // 0: car's unique ID, 
            // 1: car's x position in map coordinates, 
            // 2: car's y position in map coordinates, 
            // 3: car's x velocity in m/s, 
            // 4: car's y velocity in m/s, 
            // 5: car's s position in frenet coordinates, 
            // 6: car's d position in frenet coordinates. 

            // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
            // Later we will interpolate these waypoints with a spline and fill it in with more points that control speed.
            std::vector<double> ptsx;
            std::vector<double> ptsy;
            int prev_size = previous_path_x.size();

            bool reduce_ego_speed = false;

            double vehicle_velocity_x;
            double vehicle_velocity_y;
            double vehicle_speed;
            double vehicle_check_s;
            double vehicle_d;
            double target_vehicle_speed;
            double max_speed = MAX_SPEED;


            if(prev_size > 0)
            {
              car_s = end_path_s;
            }

            // Lane change control.
            // Need to consider speed before a lane change.

            switch(state)
            {
              case KeepState:
              {
                // Only display the warning once.
                bool warn_msg = false;
                
                // Try to increase speed just incase the vehicle in front also speeds up.
                // The sensor fusion will keep me in check if we get too close.
                speed_state = IncreaseSpeed;

                for(int i = 0; i < sensor_fusion.size(); i++)
                {
                  // Check if the vehicle is in our lane
                  vehicle_d = sensor_fusion[i][6];

                  if(CheckLaneBoundary(vehicle_d, current_lane))
                  {
                    vehicle_velocity_x = sensor_fusion[i][3];
                    vehicle_velocity_y = sensor_fusion[i][4];
                    vehicle_speed = sqrt(pow(vehicle_velocity_x,2)+pow(vehicle_velocity_y,2));
                    vehicle_check_s = sensor_fusion[i][5];

                    // Determine the vehicles position in the future
                    vehicle_check_s += ((double)prev_size*0.02*vehicle_speed);

                    if((vehicle_check_s > car_s) && 
                      ((vehicle_check_s - car_s) < 20))
                    {
                      // Approaching the car in front of us, slow down to meet the speed of the vehicle in front.
                      speed_state = ReduceSpeed;
                      if(warn_msg)
                      {
                        warn_msg = false;
                        cout << "Consider changing lanes!!!" << endl;
                      }

                      target_vehicle_speed = vehicle_speed*2;

                      state = ChangeLanesState;
                    }
                  }
                }
                break;
              }
              case ChangeLanesState:
              {
                // Remember our lane state
                int previous_lane = current_lane;

                // Flags to indicate 
                bool RightFrontClearFlag = false;
                bool RightBehindClearFlag = false;
                bool LeftFrontClearFlag = false;
                bool LeftBehindClearFlag = false;
                double value_r_min, value_r_max;
                double value_l_min, value_l_max;
                std::vector<double>::iterator result;
                std::vector<double> min_distance_l;
                std::vector<double> min_distance_r;
                std::vector<double> max_distance_l;
                std::vector<double> max_distance_r;

                for(int i = 0; i < sensor_fusion.size(); i++)
                {

                  vehicle_velocity_x = sensor_fusion[i][3];
                  vehicle_velocity_y = sensor_fusion[i][4];
                  vehicle_speed = sqrt(pow(vehicle_velocity_x,2)+pow(vehicle_velocity_y,2));
                  vehicle_check_s = sensor_fusion[i][5];
                  // Check if the vehicle is in our lane
                  vehicle_d = sensor_fusion[i][6];

                  // Determine the vehicles position in the future
                  vehicle_check_s += ((double)prev_size*0.02*vehicle_speed);

                  // Collect the distances between the ego vehicle and the other vehicle.
                  if(CheckLaneBoundary(vehicle_d, current_lane-1))
                  {
                    // Store the values of the cars in front the ego vehicle
                    if((vehicle_check_s > car_s))
                    {
                      min_distance_l.push_back(vehicle_check_s-car_s);
                    }
                    // Store the values of the cars behind the ego vehicle
                    if((vehicle_check_s < car_s))
                    {
                      max_distance_l.push_back(car_s-vehicle_check_s);
                    }
                      
                  }
                  else if(CheckLaneBoundary(vehicle_d, current_lane+1))
                  {
                    // Store the values of the cars in front the ego vehicle
                    if((vehicle_check_s > car_s))
                    {
                      min_distance_r.push_back(vehicle_check_s-car_s);
                    }

                    // Store the values of the cars behind the ego vehicle
                    if((vehicle_check_s < car_s))
                    {
                      // cout << vehicle_check_s-car_s << endl;
                      max_distance_r.push_back(car_s-vehicle_check_s);
                    }
                  }
                }

                // If there are cars in front (right lane), find the minimum distance
                if(!min_distance_r.empty())
                {
                  if((current_lane+1 <= 2))
                  {
                    result = std::min_element(std::begin(min_distance_r), std::end(min_distance_r));
                    value_r_min = min_distance_r[std::distance(std::begin(min_distance_r), result)];
                    // cout << "Right: Car is far in front: " << value_r_min << endl;

                    // Check that the closest car has a safe distance in front of the ego vehicle.
                    if(value_r_min > 20)
                    {
                      RightFrontClearFlag = true;
                    }
                  }
                }
                // There is no cars in front
                else
                {
                  if((current_lane+1 <= 2))
                  {
                    RightFrontClearFlag = true;
                  }
                }

                // If there are cars behind (right lane), find the minimum distance.
                if(!max_distance_r.empty())
                {
                  if((current_lane+1 <= 2))
                  {
                    result = std::min_element(std::begin(max_distance_r), std::end(max_distance_r));
                    value_r_max = max_distance_r[std::distance(std::begin(max_distance_r), result)];
                    // cout << "Right: Car is far behind: " << value_r_max << endl;

                    // Check that the closest car has a safe distance behind the ego vehicle.
                    if(value_r_max > 10)
                    {
                      RightBehindClearFlag = true;
                    }
                  }
                }
                // There is no cars behind
                else
                {
                  if((current_lane+1 <= 2))
                  {
                    RightBehindClearFlag = true;
                  }
                }

                // If there are cars in front (left lane), find the minimum distance
                if(!min_distance_l.empty())
                {
                  if((current_lane-1 >= 0))
                  {
                    result = std::min_element(std::begin(min_distance_l), std::end(min_distance_l));
                    value_l_min = min_distance_l[std::distance(std::begin(min_distance_l), result)];
                    // cout << "Left: Car is far in front: " << value_l << endl;

                    // Check that the closest car has a safe distance in front of the ego vehicle.
                    if(value_l_min > 20)
                    {
                      LeftFrontClearFlag = true;
                    }
                  }

                }
                else
                {
                  if((current_lane-1 >= 0))
                  {
                    LeftFrontClearFlag = true;
                  }
                }

                // If there are cars behind (left lane), find the minimum distance.
                if(!max_distance_l.empty())
                {
                  if((current_lane-1 >= 0))
                  {
                    result = std::min_element(std::begin(max_distance_l), std::end(max_distance_l));
                    value_l_max = max_distance_l[std::distance(std::begin(max_distance_l), result)];
                    // cout << "Left: Car is far behind: " << value_l << endl;

                    // Check that the closest car has a safe distance behind the ego vehicle.
                    if(value_l_max > 10)
                    {
                      LeftBehindClearFlag = true;
                    }
                  }
                }
                else
                {
                  if((current_lane-1 >= 0))
                  {
                    LeftBehindClearFlag = true;
                  }
                }

                // cout << "Right turn: " << value_r_min << " : " << -value_r_max << endl;
                // cout << "Left turn: " << value_l_min << " : " << -value_l_max << endl;
                cout << previous_lane << " : " << RightFrontClearFlag << " : " << RightBehindClearFlag << " : " << LeftFrontClearFlag << " : " << LeftBehindClearFlag << endl;

                if((RightFrontClearFlag && RightBehindClearFlag) && 
                   (LeftFrontClearFlag && LeftBehindClearFlag))
                {
                  // In this case the right and left lanes are clear, so 
                  if(min_distance_r.empty())
                  {
                    cout << "Right lane change in progress!!!" << endl;
                    current_lane += 1; 
                  }
                  else if (min_distance_l.empty())
                  {
                    cout << "Left lane change in progress!!!" << endl;
                    current_lane -= 1;
                  }
                  else if(max(value_r_min, value_l_min) == value_r_min)
                  {
                    cout << "Right lane change in progress!!!" << endl;
                    current_lane += 1; 
                  }
                  else
                  {
                    cout << "Left lane change in progress!!!" << endl;
                    current_lane -= 1;
                  }
                }
                else
                {
                  if (RightFrontClearFlag && RightBehindClearFlag)
                  {
                    // cout << "Right turn: " << value_r_min << " : " << value_r_max << endl;
                    cout << "Right lane change in progress!!!" << endl;
                    current_lane += 1; 
                                   
                  }

                  if (LeftFrontClearFlag && LeftBehindClearFlag)
                  {
                    // cout << "Left turn: " << value_l_min << " : " << value_l_max << endl;
                    cout << "Left lane change in progress!!!" << endl;
                    current_lane -= 1;
                    
                  }
                }

                if(current_lane == previous_lane)
                {
                  state = KeepState;
                }
                else
                {
                  cout << "Hold current state before considering a lane change!!!" << endl;
                  state = HoldLaneState;
                }
                break;
              }
              case HoldLaneState:
              {
                // Try to increase speed just incase the vehicle in front also speeds up.
                // The sensor fusion will keep me in check if we get too close.
                speed_state = IncreaseSpeed;

                for(int i = 0; i < sensor_fusion.size(); i++)
                {
                  // Check if the vehicle is in our lane
                  vehicle_d = sensor_fusion[i][6];

                  if(CheckLaneBoundary(vehicle_d, current_lane))
                  {
                    vehicle_velocity_x = sensor_fusion[i][3];
                    vehicle_velocity_y = sensor_fusion[i][4];
                    vehicle_speed = sqrt(pow(vehicle_velocity_x,2)+pow(vehicle_velocity_y,2));
                    vehicle_check_s = sensor_fusion[i][5];

                    // Determine the vehicles position in the future
                    vehicle_check_s += ((double)prev_size*0.02*vehicle_speed);

                    if((vehicle_check_s > car_s) && 
                      ((vehicle_check_s-car_s) < 20))
                    {
                      // Approaching the car in front of us, slow down to meet the speed of the vehicle in front.
                      speed_state = ReduceSpeed;

                      target_vehicle_speed = vehicle_speed*2;
                    }
                  }
                }

                if(loop_count < 100)
                {
                  loop_count += 1;
                }
                else
                {
                  state = KeepState;
                  loop_count = 0;
                }

                break;
              }
              default:
              {
                cout << "Error!" << endl;
              }
            }




            // Speed control
            switch(speed_state)
            {
              case IncreaseSpeed:
              {
                if(ref_vel < 49.5)
                {
                  ref_vel += 0.224;
                }
                else
                {
                  speed_state = MaintainSpeed;
                }
                break;
              }
              case ReduceSpeed:
              {

                if(target_vehicle_speed < car_speed)
                {
                  ref_vel -= 0.324;
                }
                else
                {
                  speed_state = MatchLaneSpeed;
                }
                break;
              }
              case MatchLaneSpeed:
              {
                if(car_speed < target_vehicle_speed)
                {
                  ref_vel += 0.224;
                }
                else
                {
                  speed_state = MaintainSpeed;
                }
                break;
              }
              case MaintainSpeed:
              {
                // Do nothing
                break;
              }
              default:
              {
                break;
              }
            }


            // Reference x,y, and yaw states
            // Either we will reference the starting point as where the car is or at the previous paths end point.
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            // If previous size is almost empty, use the ego vehicle as starting reference
            if(prev_size < 2)
            {
              // Use two points that make the path tangent to the car.
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
            }
            // use the previous path's end point as starting reference.
            else
            {
              //In this case, the reference state includes the previous paths end points
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double ref_x_prev = previous_path_x[prev_size-2];
              double ref_y_prev = previous_path_y[prev_size-2];
              ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

              //use two points that make the path tangent to the previous path's end point
              ptsx.push_back(ref_x_prev);
              ptsx.push_back(ref_x);

              ptsy.push_back(ref_y_prev);
              ptsy.push_back(ref_y);
            }

            // Add evenly 30m spaced points (in frenet) ahead of the starting reference.
            std::vector<double> next_wp0 = getXY(car_s+30, (2+4*current_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
            std::vector<double> next_wp1 = getXY(car_s+60, (2+4*current_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
            std::vector<double> next_wp2 = getXY(car_s+90, (2+4*current_lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);

            // Add the 3 new x points to the list
            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            // Add the 3 new y points to the list
            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            // 
            for(int i = 0; i < ptsx.size(); i++)
            {
              // Shift the cars reference angle to 0 deg.
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = (shift_x*cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw));
              ptsy[i] = (shift_x*sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw));

            }

            // Create a spline function
            tk::spline spl_function;

            // Set the x and y points to the spline.
            spl_function.set_points(ptsx, ptsy);

            // Define the next x and y points we will use for the planner.
            std::vector<double> next_x_vals;
            std::vector<double> next_y_vals;

            // Start with all of the previous path points from last time.
            for(int i = 0; i < previous_path_x.size(); i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);

            }

            // Calculate how to break up spline points so that we travel at our desired reference velocity.
            double target_x = 30.0;
            double target_y = spl_function(target_x);
            double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

            double x_add_on = 0;

            // Fill up the rest of our path planner after filling it with previous points, here we will always output 40 points
            for(int i = 1; i < 40-previous_path_x.size(); i++)
            {
              // Every 20 ms the car moves to the next point on the list
              // Number of points (N) * 20ms * velocity(mph) = distance
              // N = distance/(0.02 * velocity/2.24)
              // 2.24 is to convert mph to m/s
              double N = (target_dist/(0.02*ref_vel/2.24));
              // Divide up the points
              double x_point = x_add_on+(target_x)/N;
              // Get the y value associated with each point.
              double y_point = spl_function(x_point);

              // Offset for the next iteration 
              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              // Rotate back to normal after rotating it earlier
              x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
              y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
