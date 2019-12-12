//
// Created by yao on 19-12-10.
//
#include <iostream>
#include "ros/ros.h"
using namespace std;

int main(int argc, char **argv) {

    std::vector<int> time;

    ros::init(argc, argv, "ipopt_node");
    ros::NodeHandle nh;
    ros::Rate loop_rate(10);

    while (ros::ok()) {
        clock_t t_start = clock();

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        std::cout << "cost time: " << time.back() << std::endl;
        loop_rate.sleep();
    }
    return 0;
}
