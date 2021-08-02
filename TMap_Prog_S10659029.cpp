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

	//Ū��
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

	//��ܿ�J���I�P�g�n��
	cout << setw(10) <<"���I" << setw(10) << "�n��" << setw(10) << "�g��" << endl;
	for (int i = 0; i < p; i++)
	{
		cout << setw(12) << point[i].name << ": "<< setw(8) << point[i].latitude << setw(10) << point[i].longitude << endl;
	}
	
	//�D�X�U�I�Z���ðO��
	for (int i = 0; i < p - 1; i++)
	{
		for (int j = i + 1; j < p; j++)
		{
			distance = Distance(point[i].longitude, point[j].longitude, point[i].latitude, point[j].latitude);
			
			AdjacencyMatrix[i][j] = distance;
			AdjacencyMatrix[j][i] = distance;
		}
	}

	//���AdjacencyMatrix
	cout << endl << "AdjacencyMatrix:" << endl;
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			cout << AdjacencyMatrix[i][j] << "	";
		}
		cout << endl;
	}

	//brute-force�D�̵u�Z��
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
	//��ܵ��G
	cout << endl << "�Ĥ@�D" << endl << "�̵u���|��" << point[minPoint[0]].name << "�M" << point[minPoint[1]].name << "����, �Z��" << minDis << "����" << endl;
	
	//(2)
	//ConvexHull�t��k
	ConvexHull(AdjacencyMatrix);
	cout << endl << "�ĤG�D" << endl << "ConvexHull�򦨪����I: " << endl;
	for (int i = 0; i < ConvexHull_point.size(); i++)
	{
		cout << i + 1 << "." << ConvexHull_point[i].name << "(" << ConvexHull_point[i].latitude << ", " << ConvexHull_point[i].longitude << ")" << endl;
	}
	cout << endl;

	//(2)
	//���n
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//�p��ɶ�

	ConvexHull_Area();
	QueryPerformanceCounter(&ThatTime);
	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;

	//(2)
	//�̻��Z��
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//�p��ɶ�

	ConvexHull_Farthest();

	QueryPerformanceCounter(&ThatTime);
	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;

	//(3)
	//(a)�̵u��{
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//�p��ɶ�
	int* num;
	num = new int[p];
	for (int i = 0; i < p; i++)
	{
		num[i] = i;
	}
	ShortestTravel(0, point.size() - 1, num);
	cout << endl << "�ĤT�D" << endl;
	cout << "(a)" << endl <<"Brute Force�̵u�Z����: " << comparedistance << "����" << endl;
	QueryPerformanceCounter(&ThatTime);
	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;

	//(b).iii Convex-Hull-TSP Algorithm
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//�p��ɶ�
	ConvexHull_TSP_distance = ConvexHullTSP();
	QueryPerformanceCounter(&ThatTime);
	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;

	//(b).iiii �Z�����
	cout << endl << "(c)" << endl << "Brute Force �P ConvexHull TSP �t" << ConvexHull_TSP_distance - comparedistance << "����" << endl;

	//(4)
	//�D�X�Ҧ���{�i��
	QueryPerformanceFrequency(&PinTime);
	QueryPerformanceCounter(&ThisTime);		//�p��ɶ�
	cout << endl << "�ĥ|�D" << endl;
	Clique();
	QueryPerformanceCounter(&ThatTime);
	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;

	/*cout << endl << "���浲�G:" << endl;
	QueryPerformanceFrequency(&PinTime);

	QueryPerformanceCounter(&ThisTime);		//���

	QueryPerformanceCounter(&ThatTime);

	cout << "����ɶ��G " << (double)(ThatTime.QuadPart - ThisTime.QuadPart) * 1000 / PinTime.QuadPart << "�@��" << endl;		//��ܰ���ɶ�*/

	system("PAUSE");
	return 0;
}

int Distance(double x1, double x2, double y1, double y2) {
	double temp, width, length;

	width = abs(y1 - y2);
	length = abs(x1 - x2);
	temp = width;
	//�N���I���n�׮t���⦨�Z��
	width = int(temp * 10) % 10 * 11075 + int(temp * 100) % 10 * 1107 + int(temp * 1000) % 10 * 110 + int(temp * 10000) % 10 * 11 + int(temp * 100000) % 10;
	temp = length;
	//�N���I���g�׮t���⦨�Z��
	length = int(temp * 10) % 10 * 10175 + int(temp * 100) % 10 * 1017 + int(temp * 1000) % 10 * 101 + int(temp * 10000) % 10 * 10 + int(temp * 100000) % 10;
	
	//�β���w�z�D�X�Z��
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

			//�D�X���I�Φ����u
			m = (point[i].latitude - point[j].latitude) / (point[i].longitude - point[j].longitude); //�ײv
			constant = (-1 * m) * point[i].longitude + point[i].latitude;//�`�Ƥ����
			//�O���U�I�b�u�W�B�u������Υk��
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

			//�ˬd�I�O�_�Ҧb�P�@��
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

			//�p�G�I�Ҧb�P���A�����Φ��u�����I
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

	//����@�Ӱ_�l�I
	ConvexHull_point.push_back(temp_point[n]);

	//�N�I�Ӷ��ǩ�Jvector
	do
	{
		//�ˬd���I�s���u���t�~�@���I
		for (int i = 0; i < temp1.size(); i++)
		{
			same_point = false;
			if (ConvexHull_point[n].name == temp2[i].name)
			{
				for (int j = 0; j < ConvexHull_point.size(); j++)
				{
					//�ˬd���I�O�_�w�s�b��vector��
					if (ConvexHull_point[j].name == temp1[i].name)
						same_point = true;
				}
				if (same_point == false) //���s�b�h��J
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
					//�ˬd���I�O�_�w�s�b��vector��
					if (ConvexHull_point[j].name == temp2[i].name)
						same_point = true;
				}
				if (same_point == false) //���s�b�h��J
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
	//�Q�Φ�C���D�X�Y�߳�X���n,�I���Ӷ���
	for (int i = 0; i < ConvexHull_point.size(); i++)
	{
		//�̫�@�I�P�_�I
		if (i == (ConvexHull_point.size() - 1)) {
			sum += (ConvexHull_point[i].longitude * ConvexHull_point[0].latitude) - (ConvexHull_point[0].longitude * ConvexHull_point[i].latitude);
		}
		else{
			sum += (ConvexHull_point[i].longitude * ConvexHull_point[i + 1].latitude) - (ConvexHull_point[i + 1].longitude * ConvexHull_point[i].latitude);
		}
	}
	//��ܵ��G
	cout <<"ConvexHull �ҳ򦨳̤j���n��: " << abs(sum / 2) << "������" << endl;
}

void ConvexHull_Farthest() {
	int farthestPoint[2], temp = 0, distance;
	//�Q��brute force�D ConvexHull �����I�Z���̻�
	for (int i = 0; i < ConvexHull_point.size() - 1; i++)
	{
		for (int j = i + 1; j < ConvexHull_point.size(); j++)
		{
			distance = Distance(ConvexHull_point[i].longitude, ConvexHull_point[j].longitude, ConvexHull_point[i].latitude, ConvexHull_point[j].latitude);
			if (distance > temp) //�P�_
			{
				temp = distance;
				farthestPoint[0] = i;
				farthestPoint[1] = j;
			}
		}
	}
	//��ܵ��G
	cout << "ConvexHull �̻������I��: " << ConvexHull_point[farthestPoint[0]].name << " �M " << ConvexHull_point[farthestPoint[1]].name << ", �Z����" << temp << "����" << endl;
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
		//�I��ConvexHull�������I
		if (ConvexHull_allPoint[i] == 0)
		{
			//���I��ConvexHull���䪺�Z�����
			for (int j = 0; j < ConvexHull_point.size(); j++)
			{
				if (j != (ConvexHull_point.size() - 1))
				{
					//��X�G���@����{������
					m = ((ConvexHull_point[j].latitude - ConvexHull_point[j + 1].latitude) / (ConvexHull_point[j].longitude - ConvexHull_point[j + 1].longitude));
					b = ConvexHull_point[j].latitude - m * ConvexHull_point[j].longitude;
					//�D�I�u�W�������Z��
					PtoL = abs(m * point[i].longitude - point[i].latitude + b) / sqrt(pow(m, 2) + pow(-1, 2));
					//����I��ConvexHull��������Z���̪�ðO���U��
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
				else //�̫�@�I�P�Ĥ@�I
				{
					//��X�G���@����{������
					m = ((ConvexHull_point[j].latitude - ConvexHull_point[0].latitude) / (ConvexHull_point[j].longitude - ConvexHull_point[0].longitude));
					b = ConvexHull_point[j].latitude - m * ConvexHull_point[j].longitude;
					//�D�I�u�W�������Z��
					PtoL = (abs(m * point[i].longitude - point[i].latitude + b)) / sqrt(pow(m, 2) + pow(-1, 2));
					//����I��ConvexHull��������Z���̪�ðO���U��
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
			ConvexHull_allPoint[i] = -1;	//��ConvexHull�u�W���I�]��-1
		}

	}
	
	/*for (int i = 0; i < p; i++)
	{
		cout << ConvexHull_allPoint[i] << " ";
	}*/

	//�P�_�é�J�ﴡ�bConvexHull�������I
	do
	{
		//�̧ǥ���ConvexHull�W���@�I
		ConvexHull_TSP.push_back(ConvexHull_point[k]);

		midpoint = false;//�ΨӧP�_ConvexHull��U���I�����O�_�ݴ��J�������I
		n.clear(); //�M�żȦs���I

		for (int i = k; i < p; i++)
		{
			//�P�_�����I�bConvexHull���I��s�J�Ȧs�}�C
			if (ConvexHull_allPoint[i] == k)
			{
				place = k;
				n.push_back(point[i]);
				midpoint = true;
			}
		}

		if (midpoint)	//ConvexHull�ݬﴡ�������I
		{
			int number = n.size(); //�����ݬﴡ�X���I

			//���̶Z���Ѫ�컷�Ƨ�
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

			//�N�������I���J
			for (int i = 0; i < number; i++)
			{
				ConvexHull_TSP.push_back(n[i]);
			}
		}
		//����ConvexHull�W���I�O�_���H�O����
		k++;

	} while (k < ConvexHull_point.size() );

	//���ConvexHull TSP���G
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

	cout << endl << "ConvexHull TSP ���Z����: " << distance << "����" << endl;
	return distance;
}

void Clique() {
	int n;
	int distance = 0;
	double input_km;
	int amount = 0;
	vector<Point> G;
	vector <int> order;


	//�D�X�����I,�ín�D�ϥΪ̿�J
	/*G = polygon_centroid();
	cout << endl << "���]�����I��: ";
	cout << G.name << "(" << G.longitude << ", " << G.latitude << ")" << endl;*/
	cout << "��J�ܤֻݸg�L�X�Ӵ��I(�̤j�Ȭ�" << point.size() << "): ";
	cin >> n;
	cout << "��J�����{����W�L�h�֤���: ";
	cin >> input_km;
	input_km = input_km * 1000;
	
	//���զX C(10, n)
	do
	{
		for (int i = 0; i <= (p - n) && n > 0; ++i)
		{
			c_recur(i + 1, n - 1, p, order);
		}
		n++;
	} while (n < p + 1);

	//�զX���G
	for (int i = 0; i < path.size(); i++)
	{
		distance = 0;

		//���D���զX������G
		for (int j = 0; j < path[i].Path.size(); j++)
		{
			G.push_back(polygon_centroid(path[i].Path));
		}

		//�p��U���|�`�Z��
		for (int j = 0; j < path[i].Path.size(); j++)
		{
			if (j == 0)
				distance += Distance(G[i].longitude, point[path[i].Path[j]].longitude, G[i].latitude, point[path[i].Path[j]].latitude);
			else if(j != path[i].Path.size() - 1)
				distance += Distance(point[path[i].Path[j - 1]].longitude, point[path[i].Path[j]].longitude, point[path[i].Path[j - 1]].latitude, point[path[i].Path[j]].latitude);
			else if (j == path[i].Path.size() - 1)
				distance += Distance(G[i].longitude, point[path[i].Path[j]].longitude, G[i].latitude, point[path[i].Path[j]].latitude);		
		}

		//��ܶZ���p��ϥΪ̿�J�����G
		if (distance < input_km)
		{
			cout << endl << "�����|�����I��: ";
			cout << G[i].name << "(" << G[i].longitude << ", " << G[i].latitude << ")" << endl;
			cout << "���|��:" << endl;
			cout << G[i].name << "-->";
			for (int j = 0; j < path[i].Path.size(); j++)
			{
				cout << point[path[i].Path[j]].name << "-->";
			}
			cout << G[i].name << " ";
			cout << "�`�Z��: " << distance << "����" << endl << endl;
			amount++;
		}

	}
	cout << "�`�@��" << amount << "�Ӥ�k" << endl;
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

	G.name = "������";
	G.longitude = cx / 3 / w;
	G.latitude = cy / 3 / w;
	return G;
}

double cross(Point v1, Point v2)
{
	// �S�����k�A���q�קK�~�t�C
	return v1.longitude * v2.latitude - v1.latitude * v2.longitude;
}