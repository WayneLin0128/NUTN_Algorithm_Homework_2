#include <iostream>
#include<fstream>
#include <Windows.h>
#include <vector>
#include <string>
#include<cmath>
#include<math.h>
#include <iomanip>
#include <algorithm>
using namespace std;

struct Point {
	string name;
	double latitude;
	double longitude;
};

struct Order {
	vector <int> Path;
};

vector <Point> point;
vector<Point> ConvexHull_point;
vector<Point> ConvexHull_TSP;
vector <Order> path;

Order temp;
int p, k, r;
int comparedistance = 0;
int** AdjacencyMatrix;
int* ConvexHull_allPoint;
int Distance(double x1, double x2, double y1, double y2);
int ConvexHull(int**);
void ConvexHull_Area();
void ConvexHull_Farthest();
void ShortestTravel(int a, int b, int* num);
int ConvexHullTSP();
void Clique();
void c_recur(int k, int n, int m, vector<int> list);
Point polygon_centroid(vector<int>input_path);
double cross(Point v1, Point v2);

int main()
{
	LARGE_INTEGER ThisTime, ThatTime, PinTime;
	ifstream fin("Map2.txt");
	Point point_struct;
	int n = 0;
	int distance, minDis, minPoint[2], ConvexHull_TSP_distance;

	//讀檔
	while (!fin.eof())
	{
		fin >> point_struct.name >> point_struct.latitude >> point_struct.longitude;
		point.push_back(point_struct);
	}
	p = point.size();
	AdjacencyMatrix = new int* [p];
	for (int i = 0; i < p; i++)
	{
		AdjacencyMatrix[i] = new int[p];
	}

	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			AdjacencyMatrix[i][j] = 0;
		}
	}

	//顯示輸入的點與經緯度
	cout << setw(10) <<"景點" << setw(10) << "緯度" << setw(10) << "經度" << endl;
	for (int i = 0; i < p; i++)
	{
		cout << setw(12) << point[i].name << ": "<< setw(8) << point[i].latitude << setw(10) << point[i].longitude << endl;
	}
	
	//求出各點距離並記錄
	for (int i = 0; i < p - 1; i++)
	{
		for (int j = i + 1; j < p; j++)
		{
			distance = Distance(point[i].longitude, point[j].longitude, point[i].latitude, point[j].latitude);
			
			AdjacencyMatrix[i][j] = distance;
			AdjacencyMatrix[j][i] = distance;
		}
	}

	//顯示AdjacencyMatrix
	cout << endl << "AdjacencyMatrix:" << endl;
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			cout << AdjacencyMatrix[i][j] << "	";
		}
		cout << endl;
	}

	//brute-force求最短距離
	minDis = AdjacencyMatrix[0][1];
	for (int i = 0; i < p - 1; i++)
	{
		for (int j = i + 1; j < p; j++)
		{
			if (AdjacencyMatrix[i][j] < minDis)
			{
				minDis = AdjacencyMatrix[i][j];
				minPoint[0] = i;
				minPoint[1] = j;
			}
		}
	}

	//(1)
	//顯示結果
	cout << endl << "第一題" << endl << "最短路徑為" << point[minPoint[0]].name << "和" << point[minPoint[1]].name << "之間, 距離" << minDis << "公尺" << endl;
	
	//(2)
	//ConvexHull演算法
	ConvexHull(AdjacencyMatrix);
	cout << endl << "第二題" << endl << "ConvexHull圍成的景點: " << endl;
	for (int i = 0; i < ConvexHull_point.size(); i++)
	{
		cout << i + 1 << "." << ConvexHull_point[i].name << "(" << ConvexHull_point[i].latitude << ", " << ConvexHull_point[i].longitude << ")" << endl;
	}
	cout << endl;

	//(2)
	//面積
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//計算時間

	ConvexHull_Area();
	QueryPerformanceCounter(&ThatTime);
	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;

	//(2)
	//最遠距離
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//計算時間

	ConvexHull_Farthest();

	QueryPerformanceCounter(&ThatTime);
	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;

	//(3)
	//(a)最短行程
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//計算時間
	int* num;
	num = new int[p];
	for (int i = 0; i < p; i++)
	{
		num[i] = i;
	}
	ShortestTravel(0, point.size() - 1, num);
	cout << endl << "第三題" << endl;
	cout << "(a)" << endl <<"Brute Force最短距離為: " << comparedistance << "公尺" << endl;
	QueryPerformanceCounter(&ThatTime);
	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;

	//(b).iii Convex-Hull-TSP Algorithm
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//計算時間
	ConvexHull_TSP_distance = ConvexHullTSP();
	QueryPerformanceCounter(&ThatTime);
	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;

	//(b).iiii 距離比較
	cout << endl << "(c)" << endl << "Brute Force 與 ConvexHull TSP 差" << ConvexHull_TSP_distance - comparedistance << "公尺" << endl;

	//(4)
	//求出所有行程可能
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//計算時間
	cout << endl << "第四題" << endl;
	Clique();
	QueryPerformanceCounter(&ThatTime);
	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;

	/*cout << endl << "執行結果:" << endl;
	QueryPerformanceFrequency(&PinTime);

	QueryPerformanceCounter(&ThisTime);		//顯示

	QueryPerformanceCounter(&ThatTime);

	cout << "執行時間： " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "毫秒" << endl;		//顯示執行時間*/

	system("PAUSE");
	return 0;
}

int Distance(double x1, double x2, double y1, double y2) {
	double temp, width, length;

	width = abs(y1 - y2);
	length = abs(x1 - x2);
	temp = width;
	//將兩點的緯度差換算成距離
	width = int(temp * 10) % 10 * 11075 + int(temp * 100) % 10 * 1107 + int(temp * 1000) % 10 * 110 + int(temp * 10000) % 10 * 11 + int(temp * 100000) % 10;
	temp = length;
	//將兩點的經度差換算成距離
	length = int(temp * 10) % 10 * 10175 + int(temp * 100) % 10 * 1017 + int(temp * 1000) % 10 * 101 + int(temp * 10000) % 10 * 10 + int(temp * 100000) % 10;
	
	//用畢氏定理求出距離
	return (int(sqrt(pow(width, 2) + pow(length, 2))));
}

int ConvexHull(int** AdjacencyMatrix) {
	int* compare;
	double m, constant;
	int n = 0;
	vector<Point> temp1, temp2, temp_point;
	bool allsame = true, two_point = true, same_point = false;

	compare = new int[p];
	ConvexHull_allPoint = new int[p];

	for (int i = 0; i < p; i++)
	{
		ConvexHull_allPoint[i] = 0;
	}

	for (int i = 0; i < p - 1; i++)
	{
		for (int j = i + 1; j < p; j++) {
			allsame = true;

			//求出兩點形成的線
			m = (point[i].latitude - point[j].latitude) / (point[i].longitude - point[j].longitude); //斜率
			constant = (-1 * m) * point[i].longitude + point[i].latitude;//常數比較值
			//記錄各點在線上、線的左邊或右邊
			for (int a = 0; a < p; a++)
			{
				if (a == i || a == j)
				{
					compare[a] = 0;
				}
				else
				{
					if ((point[a].latitude - point[a].longitude * m) > constant)
						compare[a] = 1;
					else if ((point[a].latitude - point[a].longitude * m) < constant)
						compare[a] = -1;
				}

			}

			//檢查點是否皆在同一邊
			for (int b = 0; b < p; b++)
			{
				if (compare[b] > 0) {
					allsame = false;
				}
			}
			if (allsame == false) {
				allsame = true;
				for (int c = 0; c < p; c++)
				{
					if (compare[c] < 0) {
						allsame = false;
					}
				}
			}

			//如果點皆在同側，紀錄形成線之兩點
			if (allsame) {
				temp1.push_back(point[i]);
				temp2.push_back(point[j]);
				ConvexHull_allPoint[i] = 1;
				ConvexHull_allPoint[j] = 1;
				if (two_point)
				{
					temp_point.push_back(point[i]);
					two_point = false;
				}
				else
				{
					temp_point.push_back(point[j]);
				}
				
			}

		}
		two_point = true;
	}

	//先放一個起始點
	ConvexHull_point.push_back(temp_point[n]);

	//將點照順序放入vector
	do
	{
		//檢查兩點連成線的另外一個點
		for (int i = 0; i < temp1.size(); i++)
		{
			same_point = false;
			if (ConvexHull_point[n].name == temp2[i].name)
			{
				for (int j = 0; j < ConvexHull_point.size(); j++)
				{
					//檢查此點是否已存在於vector內
					if (ConvexHull_point[j].name == temp1[i].name)
						same_point = true;
				}
				if (same_point == false) //不存在則放入
				{
					ConvexHull_point.push_back(temp1[i]);
					n++;
					break;
				}
				
			}
			if (ConvexHull_point[n].name == temp1[i].name)
			{
				for (int j = 0; j < ConvexHull_point.size(); j++)
				{
					//檢查此點是否已存在於vector內
					if (ConvexHull_point[j].name == temp2[i].name)
						same_point = true;
				}
				if (same_point == false) //不存在則放入
				{
					ConvexHull_point.push_back(temp2[i]);
					n++;
					break;
				}
			}
		}

	} while (n < temp_point.size() - 1);

	return 0;
}

void ConvexHull_Area() {
	double sum = 0;
	//利用行列式求出凸殼圍出面積,點須照順序
	for (int i = 0; i < ConvexHull_point.size(); i++)
	{
		//最後一點與起點
		if (i == (ConvexHull_point.size() - 1)) {
			sum += (ConvexHull_point[i].longitude * ConvexHull_point[0].latitude) - (ConvexHull_point[0].longitude * ConvexHull_point[i].latitude);
		}
		else{
			sum += (ConvexHull_point[i].longitude * ConvexHull_point[i + 1].latitude) - (ConvexHull_point[i + 1].longitude * ConvexHull_point[i].latitude);
		}
	}
	//顯示結果
	cout <<"ConvexHull 所圍成最大面積為: " << abs(sum / 2) << "平方單位" << endl;
}

void ConvexHull_Farthest() {
	int farthestPoint[2], temp = 0, distance;
	//利用brute force求 ConvexHull 哪兩點距離最遠
	for (int i = 0; i < ConvexHull_point.size() - 1; i++)
	{
		for (int j = i + 1; j < ConvexHull_point.size(); j++)
		{
			distance = Distance(ConvexHull_point[i].longitude, ConvexHull_point[j].longitude, ConvexHull_point[i].latitude, ConvexHull_point[j].latitude);
			if (distance > temp) //判斷
			{
				temp = distance;
				farthestPoint[0] = i;
				farthestPoint[1] = j;
			}
		}
	}
	//顯示結果
	cout << "ConvexHull 最遠的兩點為: " << ConvexHull_point[farthestPoint[0]].name << " 和 " << ConvexHull_point[farthestPoint[1]].name << ", 距離為" << temp << "公尺" << endl;
}

void ShortestTravel(int a, int b, int* num) {
	int temp, firstPoint, distance = 0;
	if (a == b)
	{
		for (int i = 0; i < point.size(); i++) {
			if (i == 0) {
				temp = num[i];
				firstPoint = num[i];
			}
			
			distance += AdjacencyMatrix[temp][num[i]];
			temp = num[i];
			if (i == point.size() - 1)
			{
				distance += AdjacencyMatrix[temp][firstPoint];
				if (comparedistance == 0)
					comparedistance = distance;
				else if (comparedistance > distance)
					comparedistance = distance;
			}
		}
	}
	else
		for (int j = a; j < point.size(); j++)
		{
			swap(num[a], num[j]);
			ShortestTravel(a + 1, b, num);
			swap(num[a], num[j]);
		}
}

int ConvexHullTSP() {
	vector <Point> n;
	Point smallest;
	int k = 0, place, distance = 0;
	double m, b;
	double PtoL, temp;
	bool midpoint = false;

	for (int i = 0; i < p; i++)
	{
		//點為ConvexHull內部的點
		if (ConvexHull_allPoint[i] == 0)
		{
			//該點到ConvexHull個邊的距離比較
			for (int j = 0; j < ConvexHull_point.size(); j++)
			{
				if (j != (ConvexHull_point.size() - 1))
				{
					//算出二元一次方程式的解
					m = ((ConvexHull_point[j].latitude - ConvexHull_point[j + 1].latitude) / (ConvexHull_point[j].longitude - ConvexHull_point[j + 1].longitude));
					b = ConvexHull_point[j].latitude - m * ConvexHull_point[j].longitude;
					//求點線上的垂直距離
					PtoL = abs(m * point[i].longitude - point[i].latitude + b) / sqrt(pow(m, 2) + pow(-1, 2));
					//比較點到ConvexHull的哪個邊距離最近並記錄下來
					if (j == 0)
					{
						temp = PtoL;
						ConvexHull_allPoint[i] = j;
					}
					else if (temp > PtoL)
					{
						ConvexHull_allPoint[i] = j;
						temp = PtoL;
					}		
				}
				else //最後一點與第一點
				{
					//算出二元一次方程式的解
					m = ((ConvexHull_point[j].latitude - ConvexHull_point[0].latitude) / (ConvexHull_point[j].longitude - ConvexHull_point[0].longitude));
					b = ConvexHull_point[j].latitude - m * ConvexHull_point[j].longitude;
					//求點線上的垂直距離
					PtoL = (abs(m * point[i].longitude - point[i].latitude + b)) / sqrt(pow(m, 2) + pow(-1, 2));
					//比較點到ConvexHull的哪個邊距離最近並記錄下來
					if (temp > PtoL)
					{
						ConvexHull_allPoint[i] = j;
						temp = PtoL;
					}
				}
				
			}

		}
		else
		{
			ConvexHull_allPoint[i] = -1;	//把ConvexHull線上的點設為-1
		}

	}
	
	/*for (int i = 0; i < p; i++)
	{
		cout << ConvexHull_allPoint[i] << " ";
	}*/

	//判斷並放入穿插在ConvexHull之間的點
	do
	{
		//依序先放ConvexHull上的一點
		ConvexHull_TSP.push_back(ConvexHull_point[k]);

		midpoint = false;//用來判斷ConvexHull到下個點之間是否需插入內部的點
		n.clear(); //清空暫存的點

		for (int i = k; i < p; i++)
		{
			//判斷哪些點在ConvexHull的點後存入暫存陣列
			if (ConvexHull_allPoint[i] == k)
			{
				place = k;
				n.push_back(point[i]);
				midpoint = true;
			}
		}

		if (midpoint)	//ConvexHull需穿插內部的點
		{
			int number = n.size(); //紀錄需穿插幾個點

			//先依距離由近到遠排序
			for (int i = 0; i < number - 1; i++)
			{
				for (int j = i + 1; j < number; j++)
				{
					if (Distance(n[i].longitude, ConvexHull_point[place].longitude, n[i].latitude, ConvexHull_point[place].latitude) > Distance(n[j].longitude, ConvexHull_point[place].longitude, n[j].latitude, ConvexHull_point[place].latitude))
					{
						smallest = n[j];
						n[j] = n[i];
						n[i] = smallest;
					}
				}
			}

			//將內部的點插入
			for (int i = 0; i < number; i++)
			{
				ConvexHull_TSP.push_back(n[i]);
			}
		}
		//紀錄ConvexHull上的點是否都以記錄到
		k++;

	} while (k < ConvexHull_point.size() );

	//顯示ConvexHull TSP結果
	cout << endl << "(b)" << endl << "ConvexxHull TSP:" << endl;
	for (int i = 0; i < ConvexHull_TSP.size(); i++)
	{
		cout << ConvexHull_TSP[i].name << " ";
	}

	for (int i = 0; i < ConvexHull_TSP.size(); i++)
	{
		if (i != ConvexHull_TSP.size() - 1)
		{
			distance += Distance(ConvexHull_TSP[i].longitude, ConvexHull_TSP[i + 1].longitude, ConvexHull_TSP[i].latitude, ConvexHull_TSP[i + 1].latitude);
		}
		else
		{
			distance += Distance(ConvexHull_TSP[i].longitude, ConvexHull_TSP[0].longitude, ConvexHull_TSP[i].latitude, ConvexHull_TSP[0].latitude);
		}
	}

	cout << endl << "ConvexHull TSP 的距離為: " << distance << "公尺" << endl;
	return distance;
}

void Clique() {
	int n;
	int distance = 0;
	double input_km;
	int amount = 0;
	vector<Point> G;
	vector <int> order;


	//求出中心點,並要求使用者輸入
	/*G = polygon_centroid();
	cout << endl << "假設中心點為: ";
	cout << G.name << "(" << G.longitude << ", " << G.latitude << ")" << endl;*/
	cout << "輸入至少需經過幾個景點(最大值為" << point.size() << "): ";
	cin >> n;
	cout << "輸入此趟行程不能超過多少公里: ";
	cin >> input_km;
	input_km = input_km * 1000;
	
	//做組合 C(10, n)
	do
	{
		for (int i = 0; i <= (p - n) && n > 0; ++i)
		{
			c_recur(i + 1, n - 1, p, order);
		}
		n++;
	} while (n < p + 1);

	//組合結果
	for (int i = 0; i < path.size(); i++)
	{
		distance = 0;

		//先求此組合之重心G
		for (int j = 0; j < path[i].Path.size(); j++)
		{
			G.push_back(polygon_centroid(path[i].Path));
		}

		//計算各路徑總距離
		for (int j = 0; j < path[i].Path.size(); j++)
		{
			if (j == 0)
				distance += Distance(G[i].longitude, point[path[i].Path[j]].longitude, G[i].latitude, point[path[i].Path[j]].latitude);
			else if(j != path[i].Path.size() - 1)
				distance += Distance(point[path[i].Path[j - 1]].longitude, point[path[i].Path[j]].longitude, point[path[i].Path[j - 1]].latitude, point[path[i].Path[j]].latitude);
			else if (j == path[i].Path.size() - 1)
				distance += Distance(G[i].longitude, point[path[i].Path[j]].longitude, G[i].latitude, point[path[i].Path[j]].latitude);		
		}

		//顯示距離小於使用者輸入的結果
		if (distance < input_km)
		{
			cout << endl << "此路徑中心點為: ";
			cout << G[i].name << "(" << G[i].longitude << ", " << G[i].latitude << ")" << endl;
			cout << "路徑為:" << endl;
			cout << G[i].name << "-->";
			for (int j = 0; j < path[i].Path.size(); j++)
			{
				cout << point[path[i].Path[j]].name << "-->";
			}
			cout << G[i].name << " ";
			cout << "總距離: " << distance << "公尺" << endl << endl;
			amount++;
		}

	}
	cout << "總共有" << amount << "個方法" << endl;
}

void c_recur(int k, int n, int m, vector<int> list)
{
	list.push_back(k - 1);
	for (int i = k; i <= (m - n) && n > 0; ++i)
	{
		c_recur(i + 1, n - 1, m, list);
	}
	if (n == 0)
	{
		for (int i = 0; i < list.size(); ++i)
		{
			temp.Path.push_back(list[i]);
		}
		path.push_back(temp);
		temp.Path.clear();
	}
}

Point polygon_centroid(vector <int> input_path)
{
	Point G;
	float cx = 0, cy = 0, w = 0;
	for (int i = input_path.size() - 1, j = 0; j < input_path.size(); i = j++)
	{
		float a = cross(point[input_path[i]], point[input_path[j]]);
		cx += (point[input_path[i]].longitude + point[input_path[j]].longitude) * a;
		cy += (point[input_path[i]].latitude + point[input_path[j]].latitude) * a;
		w += a;
	}

	G.name = "停車場";
	G.longitude = cx / 3 / w;
	G.latitude = cy / 3 / w;
	return G;
}

double cross(Point v1, Point v2)
{
	// 沒有除法，儘量避免誤差。
	return v1.longitude * v2.latitude - v1.latitude * v2.longitude;
}