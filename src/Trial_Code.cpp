            // Ryans Trial

            // Sensor Fusion Parameters
            // car's unique ID, 
            // car's x position in map coordinates, 
            // car's y position in map coordinates, 
            // car's x velocity in m/s, 
            // car's y velocity in m/s, 
            // car's s position in frenet coordinates, 
            // car's d position in frenet coordinates. 

            // vector<double> next_x_vals;
            // vector<double> next_y_vals;

            // std::cout << sensor_fusion << "\n" << "\n" << std::endl;

            // std::vector<double> jmt_s, jmt_d;
            // std::vector<double> sf;
            // double previous_s = std::numeric_limits<double>::max();
            // int id;
            // double car_x_pos;
            // double car_y_pos;
            // double car_vx_pos;
            // double car_vy_pos;
            // double car_s_pos;
            // double car_d_pos;
            // double check_speed;
            // double s_diff;
            // double d_diff;

            // double result_car_s;
            // double result_car_d;
            // double result_car_speed;

            // for(int i = 0; i < sensor_fusion.size(); i++)
            // {
            //   id = sensor_fusion[i][0];
            //   car_x_pos = sensor_fusion[i][1];
            //   car_y_pos = sensor_fusion[i][2];
            //   car_vx = sensor_fusion[i][3];
            //   car_vy = sensor_fusion[i][4];
            //   car_s_pos = sensor_fusion[i][5];
            //   car_d_pos = sensor_fusion[i][6];
            //   check_speed = sqrt(vx*vx+vy*vy);
            //   s_diff = car_s_pos-car_s;
            //   d_diff = car_d_pos-car_d;

            //   // Only care about the cars in our lane.
            //   if( (car_d_pos < (2+4*lane+2)) && (car_d_pos > (2+4*lane-2)) && car_s_pos > car_s)
            //   {
                
            //     cout << id << ": " << car_s << " : " << car_d << endl;
            //     cout << id << ": " << car_s_pos << " : " << car_d_pos << endl;
            //     cout << id << ": " << s_diff << ": " << d_diff << endl;
            //     // if(s_diff < previous_s)
            //     // {
            //     //   cout << id << ": " << car_s << " : " << car_d << endl;
            //     //   cout << id << ": " << car_s_pos << " : " << car_d_pos << endl;
            //     //   cout << id << ": " << s_diff << endl;
            //     //   previous_s = s_diff;
            //     //   result_car_s = car_s_pos;
            //     //   result_car_d = car_d_pos;
            //     //   farts = true;
            //     // }
            //     cout << endl;
            //   }

            // }

            // if(farts)
            // {
            //   // cout << id << ": " << s_diff << ": " << d_diff << endl;
            //   jmt_s = JMT({car_s, 0, 0},{result_car_s, 0, 0}, 5);
            //   jmt_d = JMT({car_d, 0, 0},{result_car_d, 0, 0}, 5);
            //   // cout << pto[0] << " : " << pto[1] << " : " << pto[2] << " : " << pto[3] << " : " << pto[4] << " : " << pto[5] << endl;

            //   for(int i = 1; i < 10; i++)
            //   {
            //     double val = i*0.20;
            //     double next_s = jmt_s[0] + jmt_s[1]*pow(val,1) + jmt_s[2]*pow(val,2) + jmt_s[3]*pow(val,3) + jmt_s[4]*pow(val,4) + jmt_s[5]*pow(val,5);
            //     double next_d = jmt_d[0] + jmt_d[1]*pow(val,1) + jmt_d[2]*pow(val,2) + jmt_d[3]*pow(val,3) + jmt_d[4]*pow(val,4) + jmt_d[5]*pow(val,5);
            //     // cout << id << ": " << car_s << " : " << car_d << endl;
            //     cout << id << ": " << next_s << " : " << next_d << endl;
            //     // cout << endl;
            //     std::vector<double> xy = getXY(next_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            //     next_x_vals.push_back(xy[0]);
            //     next_y_vals.push_back(xy[1]);
            //   }
            //   cout << endl;
            //   cout << endl;
            //   farts = false;
            // }

            

            // if(farts)
            // {
            //   previous_speed = car_speed;
            //   farts = false;
            // }
            // else
            // {
            //   cout << (fabs(car_speed-previous_speed)/0.02) << endl;
            //   previous_speed = car_speed;
            // }



            // Trial 3

            // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
            // Later we will interpolate these waypoints with a spline and fill it in with more points that control speed.
            std::vector<double> ptsx;
            std::vector<double> ptsy;
            int prev_size = previous_path_x.size();


            // if(prev_size > 0)
            // {
            //   car_s = end_path_s;
            // }

            // // find ref_v to use
            // bool too_close = false;

            // for(int i = 0; i < sensor_fusion.size(); i++)
            // {
            //   // car is in my lane
            //   float d = sensor_fusion[i][6];

            //   if(d < (2+4*lane+2) && d > (2+4*lane-2))
            //   {
            //     double vx = sensor_fusion[i][3];
            //     double vy = sensor_fusion[i][4];
            //     double check_speed = sqrt(vx*vx+vy*vy);
            //     double check_car_s = sensor_fusion[i][5];

            //     // if using previous points can project s value out 
            //     check_car_s += ((double)prev_size*0.02*check_speed);
            //     // check s values greater than mine and s gap
            //     if((check_car_s > car_s) && ((check_car_s-car_s) < 30))
            //     {
            //       // do some logic here, lower ref velocity so we doint crash into the car infront of us, could also flag to try to change lanes.
            //       // ref_vel = 29.5;
            //       too_close = true;
            //     }
            //   }
            // }

            // if(too_close)
            // {
            //   ref_vel -= 0.224;
            // }
            // else if(ref_vel < 49.5)
            // {
            //   ref_vel += 0.224;
            // }



            // // reference x,y, yaw states
            // // either we will reference the starting point as where the car is or at the previous paths end point.

            // double ref_x = car_x;
            // double ref_y = car_y;
            // double ref_yaw = deg2rad(car_yaw);

            // // if previous size is almost empty, use the car as starting reference

            // if(prev_size < 2)
            // {
            //   // Use two points that make the path tangent to the car.
            //   double prev_car_x = car_x - cos(car_yaw);
            //   double prev_car_y = car_y - sin(car_yaw);

            //   ptsx.push_back(prev_car_x);
            //   ptsx.push_back(car_x);

            //   ptsy.push_back(prev_car_y);
            //   ptsy.push_back(car_y);
            // }
            // // use the previous path's end point as starting reference.
            // else
            // {
            //   // redefine reference state as previous path end point
            //   ref_x = previous_path_x[prev_size-1];
            //   ref_y = previous_path_y[prev_size-1];

            //   double ref_x_prev = previous_path_x[prev_size-2];
            //   double ref_y_prev = previous_path_y[prev_size-2];
            //   ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

            //   //use two points that make the path tangent to the previous path's end point

            //   ptsx.push_back(ref_x_prev);
            //   ptsx.push_back(ref_x);

            //   ptsy.push_back(ref_y_prev);
            //   ptsy.push_back(ref_y);
            // }

            // //In Frenet add evenly 30m spaced points ahead of the starting reference.

            // std::vector<double> next_wp0 = getXY(car_s+30, (2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
            // std::vector<double> next_wp1 = getXY(car_s+60, (2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
            // std::vector<double> next_wp2 = getXY(car_s+90, (2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);

            // ptsx.push_back(next_wp0[0]);
            // ptsx.push_back(next_wp1[0]);
            // ptsx.push_back(next_wp2[0]);

            // ptsy.push_back(next_wp0[1]);
            // ptsy.push_back(next_wp1[1]);
            // ptsy.push_back(next_wp2[1]);

            // for(int i = 0; i<ptsx.size(); i++)
            // {
            //   //shift car ref angle to 0 deg.
            //   double shift_x = ptsx[i] - ref_x;
            //   double shift_y = ptsy[i] - ref_y;

            //   ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
            //   ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));

            // }

            // //create a spline

            // tk::spline s;

            // // Set (x,y) points to the spline.

            // s.set_points(ptsx, ptsy);

            // // Define the actual (x,y) points we will use for the planner.

            // std::vector<double> next_x_vals;
            // std::vector<double> next_y_vals;

            // // Start with all of the previous path points from last time.
            // for(int i = 0; i<previous_path_x.size(); i++)
            // {
            //   next_x_vals.push_back(previous_path_x[i]);
            //   next_y_vals.push_back(previous_path_y[i]);

            // }

            // // Calculate how to break up spline points so that we travel at our desired reference velocity.
            // double target_x = 30.0;
            // double target_y = s(target_x);
            // double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

            // double x_add_on = 0;

            // // Fill up the rest of our path planner after filling it with previous points, here we will always output 50 points

            // for(int i = 1; i < 50-previous_path_x.size(); i++)
            // {
            //   // N * 0.02 * velocity = distance
            //   double N = (target_dist/(0.02*ref_vel/2.24));
            //   double x_point = x_add_on+(target_x)/N;
            //   double y_point = s(x_point);

            //   x_add_on = x_point;

            //   double x_ref = x_point;
            //   double y_ref = y_point;

            //   // rotate back to normal after rotating it earlier
            //   x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
            //   y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

            //   x_point += ref_x;
            //   y_point += ref_y;

            //   next_x_vals.push_back(x_point);
            //   next_y_vals.push_back(y_point);
            // }


            // Trial 2


            // double dist_inc = 0.35;
            // for(int i = 0; i < 50; i++)
            // {
            //   double next_s = car_s + (i+1)*dist_inc;
            //   double next_d = 6;
            //   vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            //   next_x_vals.push_back(xy[0]);
            //   next_y_vals.push_back(xy[1]); 
            //   // next_x_vals.push_back(car_x + (dist_inc*i)*cos(deg2rad(car_yaw)));
            //   // next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw))); 
            // }

            // Trial 1

            // double pos_x;
            // double pos_y;
            // double angle;
            // int path_size = previous_path_x.size();

            // for(int i = 0; i < path_size; i++)
            // {
            //     next_x_vals.push_back(previous_path_x[i]);
            //     next_y_vals.push_back(previous_path_y[i]);
            // }

            // if(path_size == 0)
            // {
            //     pos_x = car_x;
            //     pos_y = car_y;
            //     angle = deg2rad(car_yaw);
            // }
            // else
            // {
            //     pos_x = previous_path_x[path_size-1];
            //     pos_y = previous_path_y[path_size-1];

            //     double pos_x2 = previous_path_x[path_size-2];
            //     double pos_y2 = previous_path_y[path_size-2];
            //     angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
            // }

            // double dist_inc = 0.5;
            // for(int i = 0; i < 50-path_size; i++)
            // {   
            //     next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
            //     next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
            //     pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
            //     pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
            // }