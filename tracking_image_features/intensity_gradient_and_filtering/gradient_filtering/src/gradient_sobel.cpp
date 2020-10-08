#include <iostream>
#include <numeric>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;

void gradientSobel() {
  // TODO: Based on the image gradients in both x and y, compute an image
  // which contains the gradient magnitude according to the equation at the
  // beginning of this section for every pixel position. Also, apply different
  // levels of Gaussian blurring before applying the Sobel operator and compare
  // the results.
  cv::Mat img;
  img = cv::imread("../images/img1.png");

  // convert image to grayscale
  cv::Mat imgGray;
  cv::cvtColor(img, imgGray, cv::COLOR_BGR2GRAY);

  // create smooth kernel
  float gauss_data[25] = {1,  4, 7, 4,  1,  4,  16, 26, 16, 4, 7, 26, 41,
                          26, 7, 4, 16, 26, 16, 4,  1,  4,  7, 4, 1};
  cv::Mat kernel = cv::Mat(5, 5, CV_32F, gauss_data);

  for (int i=0;i<25;i++) {
    gauss_data[i] /= 255.;
  }
  // apply filter
  cv::Mat imgGray_smooth;
  cv::filter2D(imgGray, imgGray_smooth, -1, kernel, cv::Point(-1, -1), 0,
               cv::BORDER_DEFAULT);
  
  // create sobel edge kernel
  float sobel_x[9] = {-1, 0, +1, -2, 0, +2, -1, 0, +1};
  float sobel_y[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  cv::Mat kernel_x = cv::Mat(3, 3, CV_32F, sobel_x);
  cv::Mat kernel_y = cv::Mat(3, 3, CV_32F, sobel_y);
  // apply filter
  cv::Mat result_x;
  cv::filter2D(imgGray_smooth, result_x, -1, kernel_x, cv::Point(-1, -1), 0,
               cv::BORDER_DEFAULT);
  cv::Mat result_y;
  cv::filter2D(imgGray_smooth, result_y, -1, kernel_y, cv::Point(-1, -1), 0,
               cv::BORDER_DEFAULT);
  cv::Mat result;
  cv::vconcat(result_x, result_y, result);
  // show result
  string windowName = "Sobel operator (x-direction)";
  cv::namedWindow(windowName, 1); // create window
  cv::imshow(windowName, result);
  cv::waitKey(0); // wait for keyboard input before continuing
}

int main() { gradientSobel(); }
