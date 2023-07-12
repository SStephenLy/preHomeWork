#if 1
#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

class Component {               //元件
public:
	string name;
	string type;
	double value = 0;
};

class Node {                        //节点
public:
	int id = 0;
	int first_d = 0, seconde_d = 0;
	Component* firstCom = NULL;
	Component* secondCom = NULL;
};

void Convert(string inPut, string outPut) {
	string line;
	ifstream inFile(inPut);        //读文件
	ofstream outFile(outPut);      //写文件
	getline(inFile, line);
	outFile << "Title: " << line << endl;     //标题行
	
	vector<Component*> storeCom;            //保存元件的地址
	vector<Node> storeNode;                 //保存节点
	bool isCreat[10] = { 0 };
	int pos[10] = { 0 };
	int a = 0, b = 0;
	int index_com = -1; int index_node = 0;

	//读取网表文件并处理保存
	while (getline(inFile, line)) {
		
		// 忽略空白行和注释行
		if (line.empty() || line[0] == '*' || line[0] == '.') {
			continue;
		}
		else {
			//读入一行
			stringstream ss(line);
			Component* p = new Component(); 
			ss >> p->name >> a >> b >> p->value;
			
			if (p->name[0] == 'V') p->type = "VSource";
			if (p->name[0] == 'R') p->type = "Resistor";
			if (p->name[0] == 'C') p->type = "Capacitor";
			storeCom.push_back(p); index_com++;

			if (!isCreat[a]) {
				Node n;
				isCreat[a] = true;
				n.id = a;
				n.first_d = 0;
				n.firstCom = storeCom[index_com];
				storeNode.push_back(n);
				pos[a] = index_node++;
			}
			else
			{
				storeNode[pos[a]].seconde_d = 0;
				storeNode[pos[a]].secondCom = storeCom[index_com];
			}

			if (!isCreat[b]) {
				Node n;
				isCreat[b] = 1;
				n.id = b;
				n.first_d = 1;
				n.firstCom = storeCom[index_com];
				storeNode.push_back(n);
				pos[b] = index_node++;
			}
			else
			{
				storeNode[pos[b]].seconde_d = 1;
				storeNode[pos[b]].secondCom = storeCom[index_com];
			}
		}
		
	}

	//向文件输出
	outFile << "datum = 0" << "        " << "lastnode = " << storeNode.size() - 1 << endl;

	for (int i = 0; i < storeNode.size(); i++) {
		outFile << "节点" << storeNode[i].id << "        " << "所连器件数：2" << endl;

		outFile <<"       " << "编号1" << "   " << "类型：" << storeNode[i].firstCom->type
			<< " " << "连接端口：" << storeNode[i].first_d << "    " << "名称：" << storeNode[i].firstCom->name
			<< "\n" << "       " <<"value：" << storeNode[i].firstCom->value << endl;

		outFile << "       " << "编号1" << "   " << "类型：" << storeNode[i].secondCom->type
			<< " " << "连接端口：" << storeNode[i].seconde_d << "    " << "名称：" << storeNode[i].secondCom->name
			<< "\n" << "       " << "value：" << storeNode[i].secondCom->value << endl;
	}
	for (auto it = storeCom.begin(); it != storeCom.end(); it++) {
		delete* it;
	}
	inFile.close();
	outFile.close();
}


int main() {
	string inPut("./test.sp");
	string outPut("./test.txt");

	Convert(inPut, outPut);

	cout << "结束" << endl;
	return 0;
}

#endif
