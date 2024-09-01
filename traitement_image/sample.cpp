/*#################################################
#                                                 #
#                                                 #
#         NDIMBIARISOA Valdo Tsiaro Hasina        #
#           valdotsiarohasina@gmail.com           #
#                     M1 MISA                     #
#                                                 #
#                                                 #
##################################################*/



#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

using namespace std;

void run();
void showMenu();
int getChoice();
int getInteger();
cv::Mat loadImage();
void blackAndWhite(cv::Mat image);
void grayscale(cv::Mat image);
void flipVertical(cv::Mat image);
void flipHorizontal(cv::Mat image);
void flipVerticalAndHorizontal(cv::Mat image);
void saltNoise(cv::Mat image);
void colorReduce(cv::Mat image);
void remapping(cv::Mat image);
void resize(cv::Mat image);
void meanFilter(cv::Mat image);;
void gaussianFilter(cv::Mat image);
void downsampling(cv::Mat image);
void nearestNeighborStrategy(cv::Mat image);
void sobelX(cv::Mat image);
void sobelY(cv::Mat image);
void sobel(cv::Mat image);
void binarySobelLow(cv::Mat image);
void binarySobelHigh(cv::Mat image);
void cannyOperator(cv::Mat image);
void detectContoursMorphologic(cv::Mat image);
void dilated(cv::Mat image);
void eroded(cv::Mat image);
void segment(cv::Mat image);
void threshold1(int &estim, int &m1, int &m2, cv::Mat img);
void contraste(cv::Mat image);


int main() {
  run();
  return 0;
}

void run() {
  bool isRun = true;
  while (isRun) {
    showMenu();
    int choice = getChoice();
    cv::Mat image;
    switch (choice) {
      case 0:
        isRun = false;
        break;
      case 1:
        image = loadImage();
        blackAndWhite(image);
        break;
      case 2:
        image = loadImage();
        grayscale(image);
        break;
      case 3:
        image = loadImage();
        flipVertical(image);
        break;
      case 4:
        image = loadImage();
        flipHorizontal(image);
        break;
      case 5:
        image = loadImage();
        flipVerticalAndHorizontal(image);
        break;
      case 6:
        image = loadImage();
        saltNoise(image);
        break;
      case 7:
        image = loadImage();
        colorReduce(image);
        break;
      case 8:
        image = loadImage();
        remapping(image);
        break;
      case 9:
        image = loadImage();
        resize(image);
        break;
      case 10:
        image = loadImage();
        meanFilter(image);
        break;
      case 11:
        image = loadImage();
        gaussianFilter(image);
        break;
      case 12:
        image = loadImage();
        downsampling(image);
        break;
      case 13:
        image = loadImage();
        nearestNeighborStrategy(image);
        break;
      case 14:
        image = loadImage();
        sobelX(image);
        break;
      case 15:
        image = loadImage();
        sobelY(image);
        break;
      case 16:
        image = loadImage();
        sobel(image);
        break;
      case 17:
        image = loadImage();
        binarySobelLow(image);
        break;
      case 18:
        image = loadImage();
        binarySobelHigh(image);
        break;
      case 19:
        image = loadImage();
        cannyOperator(image);
        break;
      case 20:
        image = loadImage();
        detectContoursMorphologic(image);
        break;
      case 21:
        image = loadImage();
        eroded(image);
        break;
      case 22:
        image = loadImage();
        dilated(image);
        break;
      case 23:
        image = loadImage();
        segment(image);
        break;
      case 24:
        image = loadImage();
        contraste(image);
        break;
    }
  }
}

void showMenu() {
  vector<string> listMenu = {
      "quitter le programme",
      "afficher l'image en noir et blanc",
      "afficher l'image avec un nuance de gris",
      "faire une rotation vertical",
      "faire une rotation horizontal",
      "faire une rotation vertical et horizontal",
      "ajouter du bruit a l'image (noise)",
      "reduire le nombre de couleur de l'image",
      "remapper l'image (modifier le pixel de l'image)",
      "redimentionner l'image",
      "filtrage des images a l'aide du filtre median",
      "filtrage d'image a l'aide du filtre gaussien",
      "sous-echantilloner une image",
      "interpolation des valeurs de pixels",
      "detection de contours sobel x image",
      "detection de contours sobel y image",
      "detection de contours sobel image",
      "detection de contour et appliquer thresold to binary sobel image low",
      "detection de contour et appliquer thresold to binary sobel image high",
      "detection des contours de l'image avec l'operateur Canny",
      "detection des contours a l'aide du filtre morphologiques",
      "erosion de l'image",
      "dilatation de l'image",
      "segmentation d'images",
      "modifier le contraste de l'image",
  };
  for (int i = 0; i < listMenu.size(); i++) {
    cout << "\t" << i << "- " << listMenu[i] << endl;
  }
}

int getChoice() {
  int number;
  while (true) {
    cout << "Veuillez entrer votre choix (entre 0 et 24): ";
    if (std::cin >> number && number >= 0 && number <= 24) {
      break;
    } else {
      std::cout << "Saisie invalide. Assurez-vous d'entrer un nombre entre 0 et 24." << std::endl;
      std::cin.clear();
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }
  return number;
}

int getInteger() {
  int number;
  while (true) {
    if (std::cin >> number && number >= 0) {
      break;
    } else {
      std::cout << "Saisie invalide. Assurez-vous d'entrer un nombre." << std::endl;
      std::cin.clear();
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      cout << "Entrer un nombre valide: ";
    }
  }
  return number;
}

cv::Mat loadImage() {
  cv::Mat image;
  string imagePath;
  while (image.empty()) {
    cout << "Entrer le nom du fichier image: ";
    cin >> imagePath;
    image = cv::imread(imagePath);
    if (image.empty()) {
      cout << "Echec du chargement de l'image. Veuillez verifier le nom du fichier et reesayer." << endl;
    }
  }
  return image;
}

void blackAndWhite(cv::Mat image) {
  cv::Mat result = image.clone();  // Make a copy of the original image

  // Convert the image to grayscale
  cv::cvtColor(image, result, cv::COLOR_BGR2GRAY);

  cv::Mat binaryImage = cv::Mat::zeros(image.size(), CV_8UC1);

  // Iterate through the image and set pixel values to 0 or 255
  for (int i = 0; i < binaryImage.rows; i++) {
    for (int j = 0; j < binaryImage.cols; j++) {
      if (image.at<uchar>(i, j) < 128) {
        binaryImage.at<uchar>(i, j) = 0;
      } else {
        binaryImage.at<uchar>(i, j) = 255;
      }
    }
  }

  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void grayscale(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::cvtColor(result, result, cv::COLOR_BGR2GRAY);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void flipVertical(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::flip(result, result, 0);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void flipHorizontal(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::flip(result, result, 1);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void flipVerticalAndHorizontal(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::flip(result, result, -1);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void saltNoise(cv::Mat image) {
  cv::Mat result = image.clone();
  cout << "Entrer le nombre de pixels que vous voulez ecraser (n > 10000): ";
  int n = getInteger();
  // int n = 10000;
  for (int k=0; k<n; k++) {
    int i = rand()%result.cols;
    int j = rand()%result.rows;
    if (result.type() == CV_8UC1) { // gray-level image
      result.at<uchar>(j,i)= 255;
    } else if (result.type() == CV_8UC3) { // color image
      result.at<cv::Vec3b>(j,i)[0]= 255;
      result.at<cv::Vec3b>(j,i)[1]= 255;
      result.at<cv::Vec3b>(j,i)[2]= 255;
    }
    i = rand()%result.cols;
    j = rand()%result.rows;
    if (result.type() == CV_8UC1) { // gray-level image
      result.at<uchar>(j,i)= 0;
    } else if (result.type() == CV_8UC3) { // color image
      result.at<cv::Vec3b>(j,i)[0]= 0;
      result.at<cv::Vec3b>(j,i)[1]= 0;
      result.at<cv::Vec3b>(j,i)[2]= 0;
    }
  }
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void colorReduce(cv::Mat image) {
  cv::Mat result = image.clone();
  int reductionFactor = 64;
  int div = 64;
  
  int nl= result.rows; // number of lines
  int nc= result.cols * result.channels(); // total number of elements per line
  uchar div2 = div >> 1; // div2 = div/2
  for (int j=0; j<nl; j++) {
    uchar* data= result.ptr<uchar>(j);
    for (int i=0; i<nc; i++) {
      // process each pixel ---------------------
			*data++= *data/div*div + div2;
			// end of pixel processing ----------------
    }
  }
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void remapping(cv::Mat image) {
  cv::Mat result = image.clone();

  cv::Mat srcX(image.rows,image.cols,CV_32F); // x-map
	cv::Mat srcY(image.rows,image.cols,CV_32F); // y-map

	// creating the mapping
	for (int i=0; i<image.rows; i++) {
		for (int j=0; j<image.cols; j++) {
			srcX.at<float>(i,j)= j;
			srcY.at<float>(i,j)= i+3*sin(j/6.0);
			// horizontal flipping
			// srcX.at<float>(i,j)= image.cols-j-1;
			// srcY.at<float>(i,j)= i;
		}
	}
	// applying the mapping
	cv::remap(image,  // source image
		    result, // destination image
			  srcX,   // x map
			  srcY,   // y map
			  cv::INTER_LINEAR); // interpolation method

  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void resize(cv::Mat image){
  cv::Mat result = image.clone();
  int newHeight, newWidth;
  cout << "Veuillez définir la nouvelle largeur et hauteur souhaitées : " << endl;
  cout <<  "Largeur : ";
  newWidth = getInteger();
  cout <<  "Longueur : ";
  newHeight = getInteger();
  resize(image, result, cv::Size(newWidth, newHeight));
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void meanFilter(cv::Mat image) {
  cv::Mat result = image.clone();
  // cout << "Entrer la taille du filtre: ";
  // int n = getInteger();
  int n = 5;

  cv::blur(result, result, cv::Size(n, n));

  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void gaussianFilter(cv::Mat image) {
  cv::Mat result = image.clone();
  // cout << "Entrer la taille du filtre: ";
  // int n = getInteger();
  int n = 5;

  cv::GaussianBlur(image,result,
		  cv::Size(n,n), // size of the filter
		  1.5);			  // parameter controlling the shape of the Gaussian


  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void downsampling(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::GaussianBlur(result, result, cv::Size(11,11), 1.75);
	// keep only 1 of every 4 pixels
	cv::Mat reduced2(result.rows/4, result.cols/4, CV_8U);
	for (int i=0; i<reduced2.rows; i++) {
    for (int j=0; j<reduced2.cols; j++) {
      reduced2.at<uchar>(i,j)= image.at<uchar>(i*4,j*4);
    }
  }
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void nearestNeighborStrategy(cv::Mat image) {
  // cv::Mat result = image.clone();
  cv::Mat result(image.rows/4, image.cols/4, CV_8U);

	for (int i=0; i<result.rows; i++)
		for (int j=0; j<result.cols; j++)
			result.at<uchar>(i,j)= image.at<uchar>(i*4,j*4);
  cv::resize(result, result, cv::Size(), 4, 4,cv::INTER_NEAREST);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void sobelX(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::Sobel(image,  // input
		result,    // output
		CV_8U,     // image type
		1, 0,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void sobelY(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::Sobel(image,  // input
		result,    // output
		CV_8U,     // image type
		0, 1,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void sobel(cv::Mat image) {
  cv::Mat result = image.clone();
  cv::Mat sobel_x;
  cv::Mat sobel_y;
  cv::Sobel(image,  // input
		sobel_x,    // output
		CV_8U,     // image type
		1, 0,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  cv::Sobel(image,  // input
		sobel_y,    // output
		CV_8U,     // image type
		0, 1,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  // Compute norm of Sobel
	cv::Sobel(image,sobel_x,CV_16S,1,0);
	cv::Sobel(image,sobel_y,CV_16S,0,1);
	//compute the L1 norm
	result = abs(sobel_x)+abs(sobel_y);
  double sobmin, sobmax;
  cv::minMaxLoc(result,&sobmin,&sobmax);
  // Conversion to 8-bit image
  // sobelImage = -alpha*sobel + 255
  cv::Mat sobelImage;
  result.convertTo(sobelImage,CV_8U,-255./sobmax,255);
  cv::imshow("original image", image);
  cv::imshow("image result", sobelImage);
  cv::waitKey(0);
}

void binarySobelLow(cv::Mat image) {
  cv::Mat result = image.clone();

  cv::Mat sobel_x;
  cv::Mat sobel_y;
  cv::Sobel(image,  // input
		sobel_x,    // output
		CV_8U,     // image type
		1, 0,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  cv::Sobel(image,  // input
		sobel_y,    // output
		CV_8U,     // image type
		0, 1,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  // Compute norm of Sobel
	cv::Sobel(image,sobel_x,CV_16S,1,0);
	cv::Sobel(image,sobel_y,CV_16S,0,1);
	//compute the L1 norm
	result = abs(sobel_x)+abs(sobel_y);
  double sobmin, sobmax;
  cv::minMaxLoc(result,&sobmin,&sobmax);
  // Conversion to 8-bit image
  // sobelImage = -alpha*sobel + 255
  cv::Mat sobelImage;
  result.convertTo(sobelImage,CV_8U,-255./sobmax,255);


  cv::Mat sobelThresholdedLow;
	cv::threshold(sobelImage, sobelThresholdedLow, 225, 255, cv::THRESH_BINARY);

  cv::imshow("original image", image);
  cv::imshow("binary sobel low", sobelThresholdedLow);
  cv::waitKey(0);
}

void binarySobelHigh(cv::Mat image) {
  cv::Mat result = image.clone();

  cv::Mat sobel_x;
  cv::Mat sobel_y;
  cv::Sobel(image,  // input
		sobel_x,    // output
		CV_8U,     // image type
		1, 0,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  cv::Sobel(image,  // input
		sobel_y,    // output
		CV_8U,     // image type
		0, 1,      // kernel specification
		3,         // size of the square kernel 
		0.4, 128); // scale and offset
  // Compute norm of Sobel
	cv::Sobel(image,sobel_x,CV_16S,1,0);
	cv::Sobel(image,sobel_y,CV_16S,0,1);
	//compute the L1 norm
	result = abs(sobel_x)+abs(sobel_y);
  double sobmin, sobmax;
  cv::minMaxLoc(result,&sobmin,&sobmax);
  // Conversion to 8-bit image
  // sobelImage = -alpha*sobel + 255
  cv::Mat sobelImage;
  result.convertTo(sobelImage,CV_8U,-255./sobmax,255);


  cv::Mat sobelThresholdedHigh;
  cv::threshold(sobelImage, sobelThresholdedHigh, 190, 255, cv::THRESH_BINARY);

  cv::imshow("original image", image);
  cv::imshow("binary sobel high", sobelThresholdedHigh);
  cv::waitKey(0);
}

void cannyOperator(cv::Mat image) {
  cv::Mat result = image.clone();
  // Apply Canny algorithm
  cv::Canny(image, // gray-level image
  result, // output contours
  125, // low threshold
  350); // high threshold
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void detectContoursMorphologic(cv::Mat image) {
  cv::Mat result;
  cv::morphologyEx(image,result,
  cv::MORPH_GRADIENT,cv::Mat());
  // Apply threshold to obtain a binary image
  int threshold= 40;
  cv::threshold(result, result,
  threshold, 255, cv::THRESH_BINARY);
  cv::imshow("original image", image);
  cv::imshow("image result", result);
  cv::waitKey(0);
}

void eroded(cv::Mat image) {
  cv::imshow("original omage",image);

	// Erode the image with the default 3x3 structuring element (SE)
	cv::Mat eroded; // the destination image
	cv::erode(image,eroded,cv::Mat());

    // Display the eroded image
	cv::imshow("Eroded Image",eroded);
  cv::waitKey(0);
}

void dilated(cv::Mat image) {
  cv::namedWindow("original image");
	cv::imshow("original image",image);

  // Dilate the image
	cv::Mat dilated; // the destination image
	cv::dilate(image,dilated,cv::Mat());

  // Display the dilated image
	cv::namedWindow("Dilated Image");
	cv::imshow("Dilated Image",dilated);
  cv::waitKey(0);
}

void segment(cv::Mat image) {
  cv::Mat result = image.clone();

  /// Compute histogram
  int bins = 256, histo[bins];
  for (int i = 0; i < bins; i++)
    histo[i] = 0;
  for (int i = 0; i < image.rows * image.cols; i++) {
    histo[image.data[i]]++;
  }
  /// Iterated threshold calculation for bimodal image
  int thres1(200), thres2(255), m1(0), m2(0), xtop(0), vtop(0);

  /// Find top peak of histogram
  for (int i = 0; i < bins; i++)
  {
    if (histo[i] > vtop)
    {
      xtop = i;
      vtop = histo[i];
    }
  }
  thres1 = xtop;
  while (abs(thres1 - thres2) > 0)
  {
    thres2 = thres1;
    threshold1(thres1, m1, m2, image);
  }

  for (int i = 0; i < result.rows * result.cols; i++) {
    if (result.data[i] < thres1)
      result.data[i] = 0;
    else
      result.data[i] = 255;
  }

  cv::threshold(result, result, thres1, 255, 0);

  cv::imshow("original image", image);
  cv::imshow("Segmented", result);
  cv::waitKey(0);
}

void contraste(cv::Mat image) {
  cv::Mat result = image.clone();

  cout << "Veuillez définir le facteur de contraste : ";
  int n = getInteger();
  result.convertTo(result, -1, n, 0);

  cv::imshow("original image", image);
  cv::imshow("Segmented", result);
  cv::waitKey(0);
}

void threshold1(int &estim, int &m1, int &m2, cv::Mat img) {
  long t0(estim), t1(255); // assume the other pick is to the left of top pick
  long i1(0), i2(0);
  while ((t0 - t1) < 2)
  {
    m1 = 0;
    m2 = 0;
    i1 = 1;
    i2 = 1;
    for (int i = 0; i < img.rows * img.cols; i++)
    {
      if (img.data[i] < t0)
      {
        m1 += img.data[i];
        i1++;
      }
      else
      {
        m2 += img.data[i];
        i2++;
      }
    }
    m1 /= i1;
    m2 /= i2;
    t0 = t1;
    t1 = (m1 + m2) / 2;
  }
  estim = t1;
}
