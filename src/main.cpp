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

// for convenience
using json = nlohmann::json;

inline double deg2rad(double x){ 
	return x * M_PI / 180; 
}
inline double rad2deg(double x){ 
	return x * 180 / M_PI; 
}

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

double eucleadian_dist(double x1, double y1, double x2, double y2){
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int getClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y){
	// infinity
	double closestLen = 1000000; 
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++){
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = eucleadian_dist(x,y,map_x,map_y);
		if(dist < closestLen){
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = getClosestWaypoint(x,y,maps_x,maps_y);
	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];
	double heading = atan2((map_y-y),(map_x-x));
	double angle = fabs(theta-heading);

  angle = min(2*M_PI - angle, angle);

  if(angle > M_PI/4){
      closestWaypoint++;
	  if (closestWaypoint == maps_x.size()){
	      closestWaypoint = 0;
	  }
  }

  return closestWaypoint;
}

vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0){
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = eucleadian_dist(x_x,x_y,proj_x,proj_y);

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = eucleadian_dist(center_x,center_y,x_x,x_y);
	double centerToRef = eucleadian_dist(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef){
		frenet_d *= -1;
	}

	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++){
		frenet_s += eucleadian_dist(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += eucleadian_dist(0,0,proj_x,proj_y);
	return {frenet_s,frenet_d};
}

vector<double> getCartesian(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )){
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-M_PI/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

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
  	float s;  // distance along the direction of the road
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

	// Start in lane 1
	int lane = 1;

	// Reference velocity to target
	double ref_vel = 0.0; // mph

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&ref_vel](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != ""){
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry"){
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
			// The data format for each car is: [ id, x, y, vx, vy, s, d]
          	auto sensor_fusion = j[1]["sensor_fusion"];

			int prev_size = previous_path_x.size();

			if (prev_size > 0){
				car_s = end_path_s;
			}

			// Lane identifiers for other cars
			bool too_close = false;
			bool car_left = false;
			bool car_right = false;

			// Find ref_v to use, see if car is in lane
			for (int i = 0; i < sensor_fusion.size(); i++){
				// Car is in my lane
				float d = sensor_fusion[i][6];

				// Identify the lane of the car in question
				int car_lane;
				if (d >= 0 && d < 4) {
					car_lane = 0;
				} else if (d >= 4 && d < 8) {
					car_lane = 1;
				} else if (d >= 8 && d <= 12) {
					car_lane = 2;
				} else {
					continue;
				}

				// Check width of lane, in case cars are merging into our lane
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx + vy*vy);
				double check_car_s = sensor_fusion[i][5];

				// If using previous points can project an s value outwards in time
				// (What position we will be in in the future)
				// check s values greater than ours and s gap
				check_car_s += ((double)prev_size*0.02*check_speed);

				int gap = 30; // m

				// Identify whether the car is ahead, to the left, or to the right
				if (car_lane == lane){
					// Another car is ahead
					too_close |= (check_car_s > car_s) && ((check_car_s - car_s) < gap);
				} else if (car_lane - lane == 1){
					// Another car is to the right
					car_right |= ((car_s - gap) < check_car_s) && ((car_s + gap) > check_car_s);
				} else if (lane - car_lane == 1){
					// Another car is to the left
					car_left |= ((car_s - gap) < check_car_s) && ((car_s + gap) > check_car_s);
				}
			}

			// Modulate the speed to avoid collisions. Change lanes if it is safe to do so (nobody to the side)
			double acc = 0.224;
			double max_speed = 49.5;
			if (too_close){
				// A car is ahead
				// Decide to shift lanes or slow down
				if (!car_right && lane < 2){
					// No car to the right AND there is a right lane -> shift right
					lane++;
				} else if (!car_left && lane > 0){
					// No car to the left AND there is a left lane -> shift left
					lane--;
				} else{
					// Nowhere to shift -> slow down
					ref_vel -= acc;
				}
			} else{
				if (lane != 1){
					// Not in the center lane. Check if it is safe to move back
					if ((lane == 2 && !car_left) || (lane == 0 && !car_right)) {
						// Move back to the center lane
						lane = 1;
					}
				}
				
				if (ref_vel < max_speed){
					// No car ahead AND we are below the speed limit -> speed limit
					ref_vel += acc;
				}
			}

			// Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
			vector<double> ptsx;
			vector<double> ptsy;

			// Reference x, y, yaw states
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);

			// If previous size is almost empty, use the car as starting reference
			if (prev_size < 2){
				// Use two points that make the path tangent to the car
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);

				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			} else{
				// Use the previous path's endpoint as starting ref
				// Redefine reference state as previous path end point

				// Last point
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];

				// 2nd-to-last point
				double ref_x_prev = previous_path_x[prev_size-2];
				double ref_y_prev = previous_path_y[prev_size-2];
				ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

				// Use two points that make the path tangent to the path's previous endpoint
				ptsx.push_back(ref_x_prev);
				ptsx.push_back(ref_x);

				ptsy.push_back(ref_y_prev);
				ptsy.push_back(ref_y);
			}

			// Using Frenet, add 30 m evenly spaced points ahead of the starting reference
			vector<double> next_wp0 = getCartesian(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getCartesian(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getCartesian(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

			ptsx.push_back(next_wp0[0]);
			ptsx.push_back(next_wp1[0]);
			ptsx.push_back(next_wp2[0]);

			ptsy.push_back(next_wp0[1]);
			ptsy.push_back(next_wp1[1]);
			ptsy.push_back(next_wp2[1]);

			for (int i = 0; i < ptsx.size(); i++){
				// Shift car reference angle to 0 degrees
				double shift_x = ptsx[i] - ref_x;
				double shift_y = ptsy[i] - ref_y;

				ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
				ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
			}

			tk::spline s;

			// Set (x,y) points to the spline
			s.set_points(ptsx, ptsy);

			// Define the actual (x,y) points we will use for the planner
			vector<double> next_x_vals;
			vector<double> next_y_vals;

			// Start with all the previous path points from last time
			for (int i = 0; i < previous_path_x.size(); i++) {
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

			// Compute break up of spline points to travel at desired ref velocity
			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));
			double x_add_on = 0;

			// Fill up the rest of the path planner to output 50 pts
			for (int i = 1; i <= 50 - previous_path_x.size(); i++) {
				double N = (target_dist/(.02*ref_vel/2.24));
				double x_point = x_add_on + (target_x) / N;
				double y_point = s(x_point);

				x_add_on = x_point;

				double x_ref = x_point;
				double y_ref = y_point;

				// Rotate back to normal after rotating it earlier
				x_point = (x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw));
				y_point = (x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw));

				x_point += ref_x;
				y_point += ref_y;

				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}

			json msgJson;
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
