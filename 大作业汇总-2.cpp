#if 1
#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <Eigen/Dense>
#include "parser_3.h"
#include <sstream>
#include "EasyBMP.h"
using namespace std;

double stripString(char* stringIn);
void printComponents(Component* compPtr);
void printNodes(Node* nodePtr, int compFlag);
char* strComponentType(Component* compPtr);
char* ComponentTypeName(Component* compPtr);  //obtain component type name
int connectNum(Component* comPtr, Node* nodePtr); //obtain port number
bool isAccurate(double delta[], int num, double accurateValue);
void Fun(double A[][30], double x[], double b[], int n);
void convertArray(double jacMat[][30], double A[][30], double result[], double y[], int number);
Eigen::VectorXd FX(double result[], int number);
Eigen::MatrixXd JA(double jacMat[][30], int number);
void newtonRaphson(double nodeValue[], int number, double jacMat[][30], double result[], double delta[]);
void homotopy(double nodeValue[], int number, double jacMat[][30], double result[], double delta[], double ratio, const double result_0[]);
void simpleTran(double Uc[], double h, double E, double C, double R, double stop_time, int& size);


int main(int argc, char* argv[]) {
    ifstream inFile;
    ofstream outFile;
    ofstream outfile;       //add another 'outfile' to finish relevant work
    NodeHead nodeList;
    CompHead compList;
    ModelHead modelList;

    // Buffers used in parsing:
    char inName[NameLength], outName[NameLength], buf[BufLength], myOutName[NameLength],
        buf1[BufLength], buf2[BufLength], buf3[BufLength], nameBuf[NameLength],
        * bufPtr, * charPtr1, * charPtr2;
    int intBuf1, intBuf2, intBuf3, intBuf4, datum = NA, eqNum = NA, specPrintJacMNA = 0;
    double douBuf1, douBuf2, douBuf3, douBuf4;
    CompType typeBuf;
    Component* compPtr, * compPtr1, * compPtr2;
    Node* nodePtr, * nodePtr1, * nodePtr2;
    Model* modelPtr;
    TranType TtypeBuf;
    EquaType eqType = Modified;

    strcpy(inName, "NOTHING");
    strcpy(outName, "NOTHING");
    strcpy(myOutName, "NOTHING");


    // process equation types:
    if (eqNum == NA) {
        while ((eqNum != 1) && (eqNum != 2)) {
            cout << "Available Equations Types Are:" << endl
                << " <1>  Nodal" << endl
                << " <2>  Modified Nodal" << endl
                << "Please enter your choice <1, 2>:" << endl;
            cin >> buf;
            eqNum = atoi(buf);
        }
        if (eqNum == 1)
            eqType = Nodal;
        else if (eqNum == 2)
            eqType = Modified;
    }

    // process input file name:
    if (!strcmp(inName, "NOTHING")) {
        cerr << "Please enter the input Spice Netlist: <QUIT to exit>" << endl;
        cin >> inName;
        if (!strcmp(inName, "QUIT")) {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
    }
    inFile.open(inName, ios::in);
    while (!inFile) {
        cerr << inName << " is an invalid input file." << endl
            << "Please enter the input Spice Netlist: <QUIT to exit>" << endl;
        cin >> inName;
        if (!strcmp(inName, "QUIT")) {
            cerr << "Program Exited Abnormally!" << endl;
            exit(0);
        }
        inFile.open(inName, ios::in);
    }

    // process output file
    if (!strcmp(outName, "NOTHING")) {
        strcpy(outName, inName);
        strtok(outName, ".");
        strcat(outName, ".Pout");
    }
    outFile.open(outName, ios::out);
    cout << endl;


    // parsing of netlist to create linked list of models (remember to reset the fstream)
    inFile.getline(buf, BufLength);       // first line of netlist is discarded
    inFile.getline(buf, BufLength);


    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }
        strcpy(buf1, buf);
        if (!strcmp(strtok(buf1, " "), ".model")) {
            strcpy(buf2, strtok(NULL, " "));
            charPtr1 = strtok(NULL, " ");
            if (!strcmp(charPtr1, "PNP"))
                TtypeBuf = PNP;
            else if (!strcmp(charPtr1, "NPN"))
                TtypeBuf = NPN;
            else if (!strcmp(charPtr1, "NMOS"))
                TtypeBuf = NMOS;
            else if (!strcmp(charPtr1, "PMOS"))
                TtypeBuf = PMOS;

            charPtr1 = strtok(NULL, " ");
            while (charPtr1 != NULL) {
                strcpy(buf3, "");
                if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                    douBuf1 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'F') && (charPtr1[2] == '=')) {
                    douBuf2 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'B') && (charPtr1[1] == 'R') && (charPtr1[2] == '=')) {
                    douBuf3 = stripString(charPtr1);
                }
                if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[2] == '=')) {
                    douBuf4 = stripString(charPtr1);
                }
                charPtr1 = strtok(NULL, " ");
            }
            modelPtr = new Model(buf2, TtypeBuf, douBuf1, douBuf2, douBuf3, douBuf4);
            modelList.addModel(modelPtr);
        }
        inFile.getline(buf, BufLength);
    }
    inFile.close();
    inFile.open(inName, ios::in);

    char model_str[9];
    //  starting of parsing by creating linked list of components
    inFile.getline(buf, BufLength);       // first line of netlist is discarded
    inFile.getline(buf, BufLength);
    while (inFile.good()) {
        if ((buf == NULL) || (*buf == '\0')) {
            inFile.getline(buf, BufLength);
            continue;
        }

        if (isalpha(*buf)) {

            //  EDIT THIS SECTION IF NEW COMPONENTS ARE ADDED!!!
            //  we could do some rearranging in this section to catch each type in order.
            switch (*buf) {
            case 'v':
            case 'V':
                typeBuf = VSource;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'i':
            case 'I':
                cout << "I" << endl;
                typeBuf = ISource;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'q':
            case 'Q':
                typeBuf = BJT;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, NA,
                    modelList.getModel(strtok(NULL, " ")), nameBuf);
                compList.addComp(compPtr);
                break;
            case 'm':
            case 'M':
                typeBuf = MOSFET;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                intBuf3 = atoi(strtok(NULL, " "));
                intBuf4 = atoi(strtok(NULL, " "));
                compPtr = new Component(typeBuf, NA, NA, intBuf1, intBuf2, intBuf3, intBuf4,
                    modelList.getModel(strtok(NULL, " ")), nameBuf);
                compList.addComp(compPtr);
                break;
            case 'r':
            case 'R':
                typeBuf = Resistor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'd':
            case 'D':
                typeBuf = Diode;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                charPtr1 = strtok(NULL, " ");
                while (charPtr1 != NULL) {
                    if ((charPtr1[0] == 'I') && (charPtr1[1] == 'S') && (charPtr1[2] == '=')) {
                        douBuf1 = stripString(charPtr1);
                    }
                    if ((charPtr1[0] == 'T') && (charPtr1[1] == 'E') && (charPtr1[4] == '=')) {
                        douBuf2 = stripString(charPtr1);
                    }
                    charPtr1 = strtok(NULL, " ");
                }
                compPtr = new Component(typeBuf, douBuf1, douBuf2, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'c':
            case 'C':
                typeBuf = Capacitor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            case 'l':
            case 'L':
                typeBuf = Inductor;
                strcpy(nameBuf, strtok(buf, " "));
                intBuf1 = atoi(strtok(NULL, " "));
                intBuf2 = atoi(strtok(NULL, " "));
                douBuf1 = atof(strtok(NULL, " "));
                compPtr = new Component(typeBuf, douBuf1, NA, intBuf1, intBuf2, NA, NA, NULL, nameBuf);
                compList.addComp(compPtr);
                break;
            };
        }
        inFile.getline(buf, BufLength);
    }


    //  Now the components are created and it is time to set up the list of nodes.
    //  we should actually use second connector of first Source as the first Node (Datum)

    compPtr1 = compList.getComp(0);
    while (compPtr1 != NULL) {
        for (int b = 0; b < 3; b++) { /* ~> J. Erik Melo note: A component can have until 4 connectors. But here just 3 are been considered. It should change the condition to 'b <= 3' or 'b < 4'?*/
            if ((!compPtr1->isCon(b)) && (compPtr1->getConVal(b) != NA)) { //~> verify if the connector 'b' is not set && if the name of the node to which this same connector 'b' is connected is a valid name as found in the circuit file. That is, if the name is not NA, that is, if this connector was named in the instantiation of the component.
                intBuf1 = compPtr1->getConVal(b); // ~> getting the connector number as in the netlist file
                nodePtr1 = nodeList.addNode();
                nodePtr1->setNameNum(intBuf1);  // ~> naming the node as in the netlist file
                compPtr1->connect(b, nodePtr1); // ~> connecting the 'connector' of component to the node
                nodePtr1->connect(b, compPtr1); // ~> connecting the 'connection' of the node to the component

                // now search and connect all other appropriate connectors to this node.
                // error checking should be added to prevent duplicated, or skipped connectors.
                compPtr2 = compPtr1->getNext();
                while (compPtr2 != NULL) {
                    for (int c = 0; c < 3; c++) { //~> verifying which one of the others connectors (of components) are connected to the node above
                        if (compPtr2->getConVal(c) == intBuf1) { //~> if next component in the list of components has a connector with the same name (conNum) of the connector above, connect it to the same node.
                            compPtr2->connect(c, nodePtr1);
                            nodePtr1->connect(c, compPtr2);
                            break;                                    //~> As a component can only have one connector with the same name (connected in the same node), don't search the others and go out of the 'for' loop
                        }
                    }
                    compPtr2 = compPtr2->getNext();
                }
            }
        }
        compPtr1 = compPtr1->getNext();
    }

    //  At this point, we are done creating a representation of the circuit in memory
    //  now, we need to call each node to create and output its nodal equation.
    //  Each node will call the components attached for the individual contributions to the
    //  nodal equation.

      // verify that input datum is valid
    Boolean check = FALSE;
    if (datum != NA) {
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() == datum)
                check = TRUE;
            nodePtr = nodePtr->getNext();
        }
        if (check == FALSE) {
            cerr << "Datum value invalid!" << endl
                << "PROGRAM EXITED ABNORMALLY!" << endl;
            exit(0);
        }
    }

    // Loop to find lastnode
    nodePtr = nodeList.getNode(0); //~> getting the pointer to the first node, pointed by 'headNode'
    int lastnode = nodePtr->getNameNum();
    while (nodePtr != NULL) {
        lastnode = (nodePtr->getNameNum() > lastnode) ? nodePtr->getNameNum() : lastnode;
        nodePtr = nodePtr->getNext();
    }

    //  Loop to find the datum
    if (datum == NA) {
        nodePtr = nodeList.getNode(0);
        nodePtr1 = nodePtr->getNext();
        while (nodePtr1 != NULL) {
            if (nodePtr1->getCount() > nodePtr->getCount())
                nodePtr = nodePtr1;
            nodePtr1 = nodePtr1->getNext();
        }
        //datum = nodePtr->getNameNum();
        datum = 0; 
    }

    //=================================
    //~> Checking the component list
    //~> Comment this part to omit
    compPtr = compList.getComp(0);
    printComponents(compPtr);

    nodePtr = nodeList.getNode(0);
    printNodes(nodePtr, 1);

    //<~
    //==================================

      // output circuit information
    outFile << "%Parser V1.0" << endl;
    outFile << "%Input Spice Deck:  " << inName << endl;
    outFile << "%Equation Type:     ";
    if (eqType == Nodal)
        outFile << "NODAL" << endl;
    else if (eqType == Modified)
        outFile << "MODIFIED NODAL" << endl;
    outFile << "%Datum Node:        " << datum << endl;


    // create value table
    outFile << endl
        << "%*****************************************************************************" << endl;
    outFile << "%                      Component Values:" << endl;
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        compPtr->printVal(outFile);
        compPtr = compPtr->getNext();
    }
    outFile << endl
        << "%*****************************************************************************" << endl;


    // go down the nodal list and have components announce themselves
    outFile << endl << "%                      Circuit Equations: " << endl;
    nodePtr = nodeList.getNode(0);
    while (nodePtr != NULL) {
        if (nodePtr->getNameNum() != datum) {
            nodePtr->printNodal(outFile, datum, lastnode);
        }
        nodePtr = nodePtr->getNext();
    }

    //go down the component list and give equations for all sources
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        compPtr->specialPrint(outFile, datum);
        compPtr = compPtr->getNext();
    }

    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
    if (eqType != Modified) {
        compPtr = compList.getComp(0);
        while (compPtr != NULL) {
            compPtr->printSuperNode(outFile, datum, lastnode);
            compPtr = compPtr->getNext();
        }
    }


    // go down the node list and give additional MNA equations
    if (eqType == Modified) {
        nodePtr = nodeList.getNode(0);
        while (nodePtr != NULL) {
            if (nodePtr->getNameNum() != datum)
                nodePtr->printMNA(outFile, datum, lastnode);
            nodePtr = nodePtr->getNext();
        }
    }

    // print jacobians
    outFile << endl
        << "%*****************************************************************************" << endl;
    outFile << endl << "%                      Jacobians: " << endl;
    nodePtr1 = nodeList.getNode(0);
    while (nodePtr1 != NULL) {   //~> this loop handles the nodes not connected to a Vsource and those ones that are not the 'datum' node
        if (nodePtr1->getNameNum() != datum) {
            nodePtr2 = nodeList.getNode(0);
            while (nodePtr2 != NULL) {
                if (nodePtr2->getNameNum() != datum) {
                    nodePtr1->printJac(outFile, datum, nodePtr2, lastnode, eqType);
                }
                nodePtr2 = nodePtr2->getNext();
            }
        }
        nodePtr1 = nodePtr1->getNext();
    }

    // go down the component list and give equations for all sources
    compPtr = compList.getComp(0);
    while (compPtr != NULL) {
        nodePtr2 = nodeList.getNode(0);
        compPtr2 = compList.getComp(0);
        while (nodePtr2 != NULL) {
            if (nodePtr2->getNameNum() != datum) {
                compPtr->specialPrintJac(outFile, datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
            }
            nodePtr2 = nodePtr2->getNext();
        }
        specPrintJacMNA = 0;
        compPtr = compPtr->getNext();
    }




    // print the Jacobians for the additional MNA equations
    if (eqType == Modified) {
        nodePtr1 = nodeList.getNode(0);
        while (nodePtr1 != NULL) {
            if (nodePtr1->getNameNum() != datum) {
                nodePtr2 = nodeList.getNode(0);
                while (nodePtr2 != NULL) {
                    if (nodePtr2->getNameNum() != datum)
                        nodePtr1->printJacMNA(outFile, datum, nodePtr2, lastnode);
                    nodePtr2 = nodePtr2->getNext();
                }
            }
            nodePtr1 = nodePtr1->getNext();
        }
    }


    cout << endl;


    // ����ת��

    if (!strcmp(myOutName, "NOTHING")) {
        strcpy(myOutName, inName);
        strtok(myOutName, ".");
        strcat(myOutName, "out.txt");
    }
    outfile.open(myOutName, ios::out);


    nodePtr = nodeList.getNode(0);

    outfile << "datum = " << datum << "\t\t" << "lastnode = " << lastnode << endl;
    Connections* conPtr;
    while (nodePtr != NULL) {

        outfile << "�ڵ� " << nodePtr->getNameNum() << "\t\t" << "����������Ϊ��" << nodePtr->getCount() << endl;
        conPtr = nodePtr->getConList();
        while (conPtr != NULL) {
            outfile << "\t\t" << "��ţ� " << conPtr->comp->getcompNum() << "\t\t" << "���ͣ� " << ComponentTypeName(conPtr->comp) << "\t\t" << "���Ӷ˿ڣ�" << connectNum(conPtr->comp, nodePtr) << "\t\t";
            switch (conPtr->comp->getType()) {
            case VSource:
                outfile << "���ƣ�" << "VCC" << endl; break;
            default:
                outfile << "���ƣ�" << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << endl; break;
            }
            outfile << "\t\t" << "value:" << conPtr->comp->getVal() << endl;

            conPtr = conPtr->next;
        }


        nodePtr = nodePtr->getNext();
    }


    int flag = 0;
    cout << "��ѡ���ܣ�\n"
        << "<1>���KCL/KVL����\n"
        << "<2>NR����\n"
        << "<3>ͬ�����\n"
        << "<4>˲̬����" << endl;
    cin >> flag;

    if (flag != 4) {
        //���KCL/KVL����
        if (flag == 1) {
            outfile << endl
                << "%*****************************************************************************" << endl;


            // go down the nodal list and have components announce themselves
            outfile << endl << "  KCL/KVL ���� : " << endl;
            nodePtr = nodeList.getNode(0);
            while (nodePtr != NULL) {
                if (nodePtr->getNameNum() != datum) {
                    nodePtr->printNodal(outfile, datum, lastnode);
                }
                nodePtr = nodePtr->getNext();
            }

            //go down the component list and give equations for all sources
            compPtr = compList.getComp(0);
            while (compPtr != NULL) {
                compPtr->specialPrint(outfile, datum);
                compPtr = compPtr->getNext();
            }

            //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
            if (eqType != Modified) {
                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    compPtr->printSuperNode(outfile, datum, lastnode);
                    compPtr = compPtr->getNext();
                }
            }


            // go down the node list and give additional MNA equations
            if (eqType == Modified) {
                nodePtr = nodeList.getNode(0);
                while (nodePtr != NULL) {
                    if (nodePtr->getNameNum() != datum)
                        nodePtr->printMNA(outfile, datum, lastnode);
                    nodePtr = nodePtr->getNext();
                }
            }
        }
        //%****************************************************************************************************


        //����ʵ�־��������
        if (flag == 2) {
            int number = 0;
            double delta[30] = { 0 };                  //����N-R��ʽ�������ſɱȾ������ �� f(x)�Ľ��

            cout << "�����ʼ�ڵ������" << endl;
            cin >> number;
            cout << "�����ʼ�ڵ���ֵ:" << endl;
            for (int i = 0; i < number; i++) {        //   case1�����8����ֵ 0.7103 0.6725 10.0 0.7103 1.5 10.0 -0.0046 -0.0021
                cin >> nodeValue[i + 1];
            }

            nodePtr = nodeList.getNode(0);
            while (nodePtr != NULL) {
                if (nodePtr->getNameNum() != datum) {
                    nodePtr->printNodalMat(datum, lastnode, result);
                }
                nodePtr = nodePtr->getNext();
            }

            compPtr = compList.getComp(0);
            while (compPtr != NULL) {
                compPtr->specialPrintMat(datum, result);
                compPtr = compPtr->getNext();
            }


            //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
            if (eqType != Modified) {
                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    compPtr->printSuperNodeMat(datum, lastnode, result);
                    compPtr = compPtr->getNext();
                }
            }


            // go down the node list and give additional MNA equations
            if (eqType == Modified) {
                nodePtr = nodeList.getNode(0);
                while (nodePtr != NULL) {
                    if (nodePtr->getNameNum() != datum)
                        nodePtr->printMNAMat(datum, lastnode, result);
                    nodePtr = nodePtr->getNext();
                }
            }




            nodePtr1 = nodeList.getNode(0);
            while (nodePtr1 != NULL) {
                if (nodePtr1->getNameNum() != datum) {
                    nodePtr2 = nodeList.getNode(0);
                    while (nodePtr2 != NULL) {
                        if (nodePtr2->getNameNum() != datum) {
                            nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                        }
                        nodePtr2 = nodePtr2->getNext();
                    }
                }
                nodePtr1 = nodePtr1->getNext();
            }

            // go down the component list and give equations for all sources
            compPtr = compList.getComp(0);
            while (compPtr != NULL) {
                nodePtr2 = nodeList.getNode(0);
                compPtr2 = compList.getComp(0);
                while (nodePtr2 != NULL) {
                    if (nodePtr2->getNameNum() != datum) {
                        compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                    }
                    nodePtr2 = nodePtr2->getNext();
                }
                specPrintJacMNA = 0;
                compPtr = compPtr->getNext();
            }



            // print the Jacobians for the additional MNA equations
            if (eqType == Modified) {
                nodePtr1 = nodeList.getNode(0);
                while (nodePtr1 != NULL) {
                    if (nodePtr1->getNameNum() != datum) {
                        nodePtr2 = nodeList.getNode(0);
                        while (nodePtr2 != NULL) {
                            if (nodePtr2->getNameNum() != datum)
                                nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                            nodePtr2 = nodePtr2->getNext();
                        }
                    }
                    nodePtr1 = nodePtr1->getNext();
                }
            }



            int count = 1;
            double accurateValue;                //����
            int max = 10000;                      //����������


            cout << "���뾫��:" << endl;
            cin >> accurateValue;
            cout << "------------------output-------------------" << endl;


            newtonRaphson(nodeValue, number, jacMat, result, delta);   //�ȵ���һ��,�ĵ���һ�ε����ֵdelta



            while (!isAccurate(delta, number, accurateValue)) {    //���ݵõ���delta�ж��Ƿ����㾫��Ҫ��
                if (count == max) break;                          //�ﵽ�����������˳�

                for (int i = 0; i < number; i++) {
                    for (int j = 0; j < number; j++) {           //ÿ�ε���f(x)���ſɱȾ������ֵ����ͬ���ȳ�ʼ��
                        jacMat[i + 1][j + 1] = 0.0;
                    }
                    result[i + 1] = 0.0;
                    delta[i + 1] = 0.0;
                }
                count++;                                          //��¼��������



                nodePtr = nodeList.getNode(0);
                while (nodePtr != NULL) {
                    if (nodePtr->getNameNum() != datum) {
                        nodePtr->printNodalMat(datum, lastnode, result);
                    }
                    nodePtr = nodePtr->getNext();
                }

                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    compPtr->specialPrintMat(datum, result);
                    compPtr = compPtr->getNext();
                }


                //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                if (eqType != Modified) {
                    compPtr = compList.getComp(0);
                    while (compPtr != NULL) {
                        compPtr->printSuperNodeMat(datum, lastnode, result);
                        compPtr = compPtr->getNext();
                    }
                }


                // go down the node list and give additional MNA equations
                if (eqType == Modified) {
                    nodePtr = nodeList.getNode(0);
                    while (nodePtr != NULL) {
                        if (nodePtr->getNameNum() != datum)
                            nodePtr->printMNAMat(datum, lastnode, result);
                        nodePtr = nodePtr->getNext();
                    }
                }



                nodePtr1 = nodeList.getNode(0);
                while (nodePtr1 != NULL) {
                    if (nodePtr1->getNameNum() != datum) {
                        nodePtr2 = nodeList.getNode(0);
                        while (nodePtr2 != NULL) {
                            if (nodePtr2->getNameNum() != datum) {
                                nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                            }
                            nodePtr2 = nodePtr2->getNext();
                        }
                    }
                    nodePtr1 = nodePtr1->getNext();
                }

                // go down the component list and give equations for all sources
                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    nodePtr2 = nodeList.getNode(0);
                    compPtr2 = compList.getComp(0);
                    while (nodePtr2 != NULL) {
                        if (nodePtr2->getNameNum() != datum) {
                            compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                        }
                        nodePtr2 = nodePtr2->getNext();
                    }
                    specPrintJacMNA = 0;
                    compPtr = compPtr->getNext();
                }


                // print the Jacobians for the additional MNA equations
                if (eqType == Modified) {
                    nodePtr1 = nodeList.getNode(0);
                    while (nodePtr1 != NULL) {
                        if (nodePtr1->getNameNum() != datum) {
                            nodePtr2 = nodeList.getNode(0);
                            while (nodePtr2 != NULL) {
                                if (nodePtr2->getNameNum() != datum)
                                    nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                                nodePtr2 = nodePtr2->getNext();
                            }
                        }
                        nodePtr1 = nodePtr1->getNext();
                    }
                }

                newtonRaphson(nodeValue, number, jacMat, result, delta);        //��ѭ�����������

            }

            //�˳�ѭ��������

            cout << "��������:" << "  " << count << endl;       //�����������

            cout << "���:" << endl;
            for (int i = 0; i < number; i++) {
                cout << "x(" << i + 1 << ") =    " << nodeValue[i + 1] << endl;
            }

            /////~~~~~~
            cout << "result���:" << endl;
            for (int i = 0; i < number; i++) {
                cout << "x(" << i + 1 << ") =    " << result[i + 1] << endl;
            }

        }

        if(flag==3){
            //%********************************************************************************************
            //ͬ�׷����

            double N = 0;        //[0,1]����ֶθ���
            /*number �� delta����ʹ��֮ǰ�����*/
            double result_0[30] = { 0 };
            double ratio = 0;
            int number = 0;
            double delta[30] = { 0 };
            double Gleak = 1e-3;
            double a[8] = { 0.5000, 0.4799, 0.9047, 0.6099, 0.6177, 0.8594, 0.8055, 0.5767 };
            double path1, lambda = 0;

            cout << "\n" << "--------------------ͬ�׷�---------------------" << endl;
            /*cout << "���벽����" << endl;
            cin >> path1;*/
            //cout << "����������ֶ���N��" << endl;
            //cin >> N;
            cout << "�����ʼ���ݸ�����" << endl;
            cin >> number;
            for (int i = 0; i < number; i++) {
                
                nodeValue[i + 1] = a[i];
                
            }

            //cout << "�����ʼ��ֵ:" << endl;
                   //   case1�����8����ֵ 0.7103 0.6725 10.0 0.7103 1.5 10.0 -0.0046 -0.0021
            for (lambda = 0; lambda < 1; lambda += 0.01) {
                for (int i = 1; i <= number; i++) {
                    result[i] = 0;

                }

                nodePtr = nodeList.getNode(0);
                while (nodePtr != NULL) {
                    if (nodePtr->getNameNum() != datum) {
                        nodePtr->printNodalMat(datum, lastnode, result);
                    }
                    nodePtr = nodePtr->getNext();
                }

                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    compPtr->specialPrintMat(datum, result);
                    compPtr = compPtr->getNext();
                }

                //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                if (eqType != Modified) {
                    compPtr = compList.getComp(0);
                    while (compPtr != NULL) {
                        compPtr->printSuperNodeMat(datum, lastnode, result);
                        compPtr = compPtr->getNext();
                    }
                }

                // go down the node list and give additional MNA equations
                if (eqType == Modified) {
                    nodePtr = nodeList.getNode(0);
                    while (nodePtr != NULL) {
                        if (nodePtr->getNameNum() != datum)
                            nodePtr->printMNAMat(datum, lastnode, result);
                        nodePtr = nodePtr->getNext();
                    }
                }

                nodePtr1 = nodeList.getNode(0);
                while (nodePtr1 != NULL) {
                    if (nodePtr1->getNameNum() != datum) {
                        nodePtr2 = nodeList.getNode(0);
                        while (nodePtr2 != NULL) {
                            if (nodePtr2->getNameNum() != datum) {
                                nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                            }
                            nodePtr2 = nodePtr2->getNext();
                        }
                    }
                    nodePtr1 = nodePtr1->getNext();
                }

                // go down the component list and give equations for all sources
                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    nodePtr2 = nodeList.getNode(0);
                    compPtr2 = compList.getComp(0);
                    while (nodePtr2 != NULL) {
                        if (nodePtr2->getNameNum() != datum) {
                            compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                        }
                        nodePtr2 = nodePtr2->getNext();
                    }
                    specPrintJacMNA = 0;
                    compPtr = compPtr->getNext();
                }

                // print the Jacobians for the additional MNA equations
                if (eqType == Modified) {
                    nodePtr1 = nodeList.getNode(0);
                    while (nodePtr1 != NULL) {
                        if (nodePtr1->getNameNum() != datum) {
                            nodePtr2 = nodeList.getNode(0);
                            while (nodePtr2 != NULL) {
                                if (nodePtr2->getNameNum() != datum)
                                    nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                                nodePtr2 = nodePtr2->getNext();
                            }
                        }
                        nodePtr1 = nodePtr1->getNext();
                    }
                }
                Eigen::MatrixXd J = JA(jacMat, number).inverse();
                J =J*lambda;
                for (int i = 0; i < number; i++) {
                    J(i, i) += (1 - lambda) * Gleak;
                }
                J = J.inverse();
                double Hk[30] = { 0 };
                for (int i = 0; i < number; i++) {
                    Hk[i] = result[i + 1]*lambda+(1-lambda)*Gleak*(nodeValue[i+1]-a[i]);
                }
                double temp[30] = { 0 };
                for (int i = 0; i < number; i++) {
                    for (int j = 0; j < number; j++) {
                        temp[i] += J(i, j) * Hk[j];
                    }
                }
                for (int i = 0; i < number; i++) {
                    nodeValue[i + 1] -= temp[i];
                }


                ////result_0����f(x0)
                //for (int i = 0; i < number; i++) {
                //    result_0[i + 1] = result[i + 1];
                //}

                ////ͬ�׵���
                //for (int k = 0; k < (int)N; k++) {

                //    ratio = k / N - 1;

                //    for (int i = 0; i < number; i++) {
                //        for (int j = 0; j < number; j++) {           //ÿ�ε���f(x)���ſɱȾ������ֵ����ͬ���ȳ�ʼ��
                //            jacMat[i + 1][j + 1] = 0.0;
                //        }
                //        result[i + 1] = 0.0;
                //        delta[i + 1] = 0.0;
                //    }

                //    //f(x)���ſɱȾ���ֵ
                //    nodePtr = nodeList.getNode(0);
                //    while (nodePtr != NULL) {
                //        if (nodePtr->getNameNum() != datum) {
                //            nodePtr->printNodalMat(datum, lastnode, result);
                //        }
                //        nodePtr = nodePtr->getNext();
                //    }

                //    compPtr = compList.getComp(0);
                //    while (compPtr != NULL) {
                //        compPtr->specialPrintMat(datum, result);
                //        compPtr = compPtr->getNext();
                //    }

                //    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                //    if (eqType != Modified) {
                //        compPtr = compList.getComp(0);
                //        while (compPtr != NULL) {
                //            compPtr->printSuperNodeMat(datum, lastnode, result);
                //            compPtr = compPtr->getNext();
                //        }
                //    }

                //    // go down the node list and give additional MNA equations
                //    if (eqType == Modified) {
                //        nodePtr = nodeList.getNode(0);
                //        while (nodePtr != NULL) {
                //            if (nodePtr->getNameNum() != datum)
                //                nodePtr->printMNAMat(datum, lastnode, result);
                //            nodePtr = nodePtr->getNext();
                //        }
                //    }

                //    nodePtr1 = nodeList.getNode(0);
                //    while (nodePtr1 != NULL) {
                //        if (nodePtr1->getNameNum() != datum) {
                //            nodePtr2 = nodeList.getNode(0);
                //            while (nodePtr2 != NULL) {
                //                if (nodePtr2->getNameNum() != datum) {
                //                    nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                //                }
                //                nodePtr2 = nodePtr2->getNext();
                //            }
                //        }
                //        nodePtr1 = nodePtr1->getNext();
                //    }

                //    // go down the component list and give equations for all sources
                //    compPtr = compList.getComp(0);
                //    while (compPtr != NULL) {
                //        nodePtr2 = nodeList.getNode(0);
                //        compPtr2 = compList.getComp(0);
                //        while (nodePtr2 != NULL) {
                //            if (nodePtr2->getNameNum() != datum) {
                //                compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                //            }
                //            nodePtr2 = nodePtr2->getNext();
                //        }
                //        specPrintJacMNA = 0;
                //        compPtr = compPtr->getNext();
                //    }

                //    // print the Jacobians for the additional MNA equations
                //    if (eqType == Modified) {
                //        nodePtr1 = nodeList.getNode(0);
                //        while (nodePtr1 != NULL) {
                //            if (nodePtr1->getNameNum() != datum) {
                //                nodePtr2 = nodeList.getNode(0);
                //                while (nodePtr2 != NULL) {
                //                    if (nodePtr2->getNameNum() != datum)
                //                        nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                //                    nodePtr2 = nodePtr2->getNext();
                //                }
                //            }
                //            nodePtr1 = nodePtr1->getNext();
                //        }
                //    }

                //    homotopy(nodeValue, number, jacMat, result, delta, ratio, result_0);

                //}

                //for (int i = 0; i < number; i++) {
                //    delta[i + 1] = 0;
                //}
                //cout << "\n" << "ͬ�׵�����ĳ�ʼֵ��" << endl;
                //for (int i = 0; i < number; i++) {
                //    cout << nodeValue[i + 1] << " ";
                //}
                //cout << endl;

                //int count = 1;
                //double accurateValue;
                //cout << "���뾫�ȣ�" << endl;
                //cin >> accurateValue;
                //int max = 4000;

                //cout << "--------------------------------output------------------------------------" << endl;


                //newtonRaphson(nodeValue, number, jacMat, result, delta);   //�ȵ���һ��,�ĵ���һ�ε����ֵdelta



                //while (!isAccurate(delta, number, accurateValue)) {    //���ݵõ���delta�ж��Ƿ����㾫��Ҫ��
                //    if (count == max) break;                          //�ﵽ�����������˳�

                //    for (int i = 0; i < number; i++) {
                //        for (int j = 0; j < number; j++) {           //ÿ�ε���f(x)���ſɱȾ������ֵ����ͬ���ȳ�ʼ��
                //            jacMat[i + 1][j + 1] = 0.0;
                //        }
                //        result[i + 1] = 0.0;
                //        delta[i + 1] = 0.0;
                //    }
                //    count++;                                          //��¼��������

                //    nodePtr = nodeList.getNode(0);
                //    while (nodePtr != NULL) {
                //        if (nodePtr->getNameNum() != datum) {
                //            nodePtr->printNodalMat(datum, lastnode, result);
                //        }
                //        nodePtr = nodePtr->getNext();
                //    }

                //    compPtr = compList.getComp(0);
                //    while (compPtr != NULL) {
                //        compPtr->specialPrintMat(datum, result);
                //        compPtr = compPtr->getNext();
                //    }


                //    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                //    if (eqType != Modified) {
                //        compPtr = compList.getComp(0);
                //        while (compPtr != NULL) {
                //            compPtr->printSuperNodeMat(datum, lastnode, result);
                //            compPtr = compPtr->getNext();
                //        }
                //    }


                //    // go down the node list and give additional MNA equations
                //    if (eqType == Modified) {
                //        nodePtr = nodeList.getNode(0);
                //        while (nodePtr != NULL) {
                //            if (nodePtr->getNameNum() != datum)
                //                nodePtr->printMNAMat(datum, lastnode, result);
                //            nodePtr = nodePtr->getNext();
                //        }
                //    }



                //    nodePtr1 = nodeList.getNode(0);
                //    while (nodePtr1 != NULL) {
                //        if (nodePtr1->getNameNum() != datum) {
                //            nodePtr2 = nodeList.getNode(0);
                //            while (nodePtr2 != NULL) {
                //                if (nodePtr2->getNameNum() != datum) {
                //                    nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                //                }
                //                nodePtr2 = nodePtr2->getNext();
                //            }
                //        }
                //        nodePtr1 = nodePtr1->getNext();
                //    }

                //    // go down the component list and give equations for all sources
                //    compPtr = compList.getComp(0);
                //    while (compPtr != NULL) {
                //        nodePtr2 = nodeList.getNode(0);
                //        compPtr2 = compList.getComp(0);
                //        while (nodePtr2 != NULL) {
                //            if (nodePtr2->getNameNum() != datum) {
                //                compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                //            }
                //            nodePtr2 = nodePtr2->getNext();
                //        }
                //        specPrintJacMNA = 0;
                //        compPtr = compPtr->getNext();
                //    }


                //    // print the Jacobians for the additional MNA equations
                //    if (eqType == Modified) {
                //        nodePtr1 = nodeList.getNode(0);
                //        while (nodePtr1 != NULL) {
                //            if (nodePtr1->getNameNum() != datum) {
                //                nodePtr2 = nodeList.getNode(0);
                //                while (nodePtr2 != NULL) {
                //                    if (nodePtr2->getNameNum() != datum)
                //                        nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                //                    nodePtr2 = nodePtr2->getNext();
                //                }
                //            }
                //            nodePtr1 = nodePtr1->getNext();
                //        }
                //    }

                //    newtonRaphson(nodeValue, number, jacMat, result, delta);        //��ѭ�����������

                //}

                //�˳�ѭ��������

                //cout << "iteration number:" << "  " << count << endl;       //�����������

                cout << "���:" << endl;
                for (int i = 0; i < number; i++) {
                    cout << "x(" << i + 1 << ") =    " << nodeValue[i + 1] << endl;
                }

            }
        }
        
    }

    //%********************************************************************************************
    
    // ˲̬����
    if (flag == 4) {
        //~~~~~~~~~~~
        //cout << "*******************************˲̬����*******************************" << endl;
        //double h = 0;                               //���沽��
        //double stop_time = 0;                       //�����ֹʱ��
        //cout << "�����벽����" << endl;
        //cin >> h;
        //cout << "���������ʱ�䣺" << endl;
        //cin >> stop_time;
        //double E = 0, C = 0, R = 0;
        //double array[100] = { 0 };
        //int size = 0;


        ////�õ�����
        //Component* comp_ptr;
        //comp_ptr = compList.getComp(0);
        //while (comp_ptr != NULL) {
        //    switch (comp_ptr->getType()) {
        //    case VSource:
        //        E = comp_ptr->getVal(); break;
        //    case Resistor:
        //        R = comp_ptr->getVal(); break;
        //    case Capacitor:
        //        C = comp_ptr->getVal(); break;
        //    }

        //    comp_ptr = comp_ptr->getNext();
        //}

        //
        //cout << endl;
        ////����ÿһ���ĵ�ѹ������array��
        //simpleTran(array, h, E, C, R, stop_time, size);

        //// ����һ��λͼ����
        //BMP image;

        //// ����ͼ���С����ɫ���
        //int width = 800;
        //int height = 600;
        //image.SetSize(width, height);
        //image.SetBitDepth(24);


        //// ���������е����ֵ����Сֵ
        //int minValue = array[0];
        //int maxValue = array[0];
        //for (int i = 1; i < width; ++i) {
        //    if (array[i] < minValue) {
        //        minValue = array[i];
        //    }
        //    if (array[i] > maxValue) {
        //        maxValue = array[i];
        //    }
        //}

        //// ѭ������ͼ���ÿ�����ز���������ɫ
        //for (int x = 0; x < width; ++x) {
        //    // ����y������
        //    int y = static_cast<int>((static_cast<double>(array[x] - minValue) / static_cast<double>(maxValue - minValue)) * height);

        //    // �������ص����ɫ
        //    image(x, height - y - 1)->Red = 255;   // ���ú�ɫ����
        //    image(x, height - y - 1)->Green = 0;   // ������ɫ����
        //    image(x, height - y - 1)->Blue = 0;    // ������ɫ����
        //}

        //// ����ͼ��ΪBMP�ļ�
        //image.WriteToFile("output.bmp");

        //cout << "ͼ���ѱ���Ϊoutput.bmp�ļ�" << endl;
        //~~~~~~~~~~~

        //���޸�
        cout << "*******************************˲̬����*******************************" << endl;
        
        int Capacitor_count = 0;
        Component* comPtr1 = compList.getComp(0);
        while (comPtr1){
            if (comPtr1->getType() == Capacitor) {
                vector<double> vs(50);
                Vs.push_back(vs);
                Vs[Capacitor_count].push_back(0);
                
                Ca[0][Capacitor_count] = comPtr1->con0.node->getNameNum();
                Ca[1][Capacitor_count] = comPtr1->con1.node->getNameNum();

                Capacitor_count++;
                
            }
            comPtr1 = comPtr1->getNext();
        }

        double stop_time = 0;                       //�����ֹʱ��
        cout << "�����벽����" << endl;
        cin >> h;
        cout << "���������ʱ�䣺" << endl;
        cin >> stop_time;
        cout << "���뾫��:" << endl;
        double accurateValue;
        cin >> accurateValue;
        cout << "�����ʼ�ڵ������" << endl;
        int number1 = 0;
        cin >> number1;
        cout << "�����ʼ�ڵ���ֵ:" << endl;
        for (int i = 0; i < number1; i++) {        //   case1�����8����ֵ 0.7103 0.6725 10.0 0.7103 1.5 10.0 -0.0046 -0.0021
            cin >> nodeValue[i + 1];
        }
        int h_count = stop_time / h;
        for (int i = 0; i < h_count; i++) {
            for (int j = 0; j < Capacitor_count; j++) {
                //
                Uk[j] = Vs[j][i];                            //������ͷ�ļ�U(k)
            }
                int number = number1;
                double delta[30] = { 0 };                  //����N-R��ʽ�������ſɱȾ������ �� f(x)�Ľ��

            
            //cout << "�����ʼ�ڵ���ֵ:" << endl;
            //for (int i = 0; i < number; i++) {        //   case1�����8����ֵ 0.7103 0.6725 10.0 0.7103 1.5 10.0 -0.0046 -0.0021
            //    cin >> nodeValue[i + 1];              //      0.71  0.67  10   0.71  1.5  10  -0.004  -0.002
            //}                                          //     0.7   0.6  10  0.7 1.5 10 

                nodePtr = nodeList.getNode(0);
                while (nodePtr != NULL) {
                    if (nodePtr->getNameNum() != datum) {
                        nodePtr->printNodalMat(datum, lastnode, result);
                    }
                    nodePtr = nodePtr->getNext();
                }

                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    compPtr->specialPrintMat(datum, result);
                    compPtr = compPtr->getNext();
                }


            //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                if (eqType != Modified) {
                    compPtr = compList.getComp(0);
                    while (compPtr != NULL) {
                        compPtr->printSuperNodeMat(datum, lastnode, result);
                        compPtr = compPtr->getNext();
                    }
                }


                // go down the node list and give additional MNA equations
                if (eqType == Modified) {
                    nodePtr = nodeList.getNode(0);
                    while (nodePtr != NULL) {
                        if (nodePtr->getNameNum() != datum)
                            nodePtr->printMNAMat(datum, lastnode, result);
                        nodePtr = nodePtr->getNext();
                    }
                }




                nodePtr1 = nodeList.getNode(0);
                while (nodePtr1 != NULL) {
                    if (nodePtr1->getNameNum() != datum) {
                        nodePtr2 = nodeList.getNode(0);
                        while (nodePtr2 != NULL) {
                            if (nodePtr2->getNameNum() != datum) {
                                nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                            }
                            nodePtr2 = nodePtr2->getNext();
                        }
                    }
                    nodePtr1 = nodePtr1->getNext();
                }

                // go down the component list and give equations for all sources
                compPtr = compList.getComp(0);
                while (compPtr != NULL) {
                    nodePtr2 = nodeList.getNode(0);
                    compPtr2 = compList.getComp(0);
                    while (nodePtr2 != NULL) {
                        if (nodePtr2->getNameNum() != datum) {
                            compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                        }
                        nodePtr2 = nodePtr2->getNext();
                    }
                    specPrintJacMNA = 0;
                    compPtr = compPtr->getNext();
                }



                // print the Jacobians for the additional MNA equations
                if (eqType == Modified) {
                    nodePtr1 = nodeList.getNode(0);
                    while (nodePtr1 != NULL) {
                        if (nodePtr1->getNameNum() != datum) {
                            nodePtr2 = nodeList.getNode(0);
                            while (nodePtr2 != NULL) {
                                if (nodePtr2->getNameNum() != datum)
                                    nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                                nodePtr2 = nodePtr2->getNext();
                            }
                        }
                        nodePtr1 = nodePtr1->getNext();
                    }
                }



                int count = 1;
                               //����
                int max = 10000;                      //����������


            
                //cout << "------------------output-------------------" << endl;


                newtonRaphson(nodeValue, number, jacMat, result, delta);   //�ȵ���һ��,�ĵ���һ�ε����ֵdelta



                while (!isAccurate(delta, number, accurateValue)) {    //���ݵõ���delta�ж��Ƿ����㾫��Ҫ��
                    if (count == max) break;                          //�ﵽ�����������˳�

                    for (int i = 0; i < number; i++) {
                        for (int j = 0; j < number; j++) {           //ÿ�ε���f(x)���ſɱȾ������ֵ����ͬ���ȳ�ʼ��
                            jacMat[i + 1][j + 1] = 0.0;
                        }
                        result[i + 1] = 0.0;
                        delta[i + 1] = 0.0;
                    }
                    count++;                                          //��¼��������



                    nodePtr = nodeList.getNode(0);
                    while (nodePtr != NULL) {
                        if (nodePtr->getNameNum() != datum) {
                            nodePtr->printNodalMat(datum, lastnode, result);
                        }
                        nodePtr = nodePtr->getNext();
                    }

                    compPtr = compList.getComp(0);
                    while (compPtr != NULL) {
                        compPtr->specialPrintMat(datum, result);
                        compPtr = compPtr->getNext();
                    }


                    //~> go down the component list and give supernode equations for all float sources (Nodal Analysis)
                    if (eqType != Modified) {
                        compPtr = compList.getComp(0);
                        while (compPtr != NULL) {
                            compPtr->printSuperNodeMat(datum, lastnode, result);
                            compPtr = compPtr->getNext();
                        }
                    }


                    // go down the node list and give additional MNA equations
                    if (eqType == Modified) {
                        nodePtr = nodeList.getNode(0);
                        while (nodePtr != NULL) {
                            if (nodePtr->getNameNum() != datum)
                                nodePtr->printMNAMat(datum, lastnode, result);
                            nodePtr = nodePtr->getNext();
                        }
                    }



                    nodePtr1 = nodeList.getNode(0);
                    while (nodePtr1 != NULL) {
                        if (nodePtr1->getNameNum() != datum) {
                            nodePtr2 = nodeList.getNode(0);
                            while (nodePtr2 != NULL) {
                                if (nodePtr2->getNameNum() != datum) {
                                    nodePtr1->printJacMat(datum, nodePtr2, lastnode, eqType, jacMat);
                                }
                                nodePtr2 = nodePtr2->getNext();
                            }
                        }
                        nodePtr1 = nodePtr1->getNext();
                    }

                    // go down the component list and give equations for all sources
                    compPtr = compList.getComp(0);
                    while (compPtr != NULL) {
                        nodePtr2 = nodeList.getNode(0);
                        compPtr2 = compList.getComp(0);
                        while (nodePtr2 != NULL) {
                            if (nodePtr2->getNameNum() != datum) {
                                compPtr->specialPrintJacMat(datum, nodePtr2, lastnode, eqType, compPtr2, &specPrintJacMNA, jacMat); // ~> specPrintJacMNA is used to verify if the jacobians w.r.t. the Modified equations was already printed to print only once.
                            }
                            nodePtr2 = nodePtr2->getNext();
                        }
                        specPrintJacMNA = 0;
                        compPtr = compPtr->getNext();
                    }


                    // print the Jacobians for the additional MNA equations
                    if (eqType == Modified) {
                        nodePtr1 = nodeList.getNode(0);
                        while (nodePtr1 != NULL) {
                            if (nodePtr1->getNameNum() != datum) {
                                nodePtr2 = nodeList.getNode(0);
                                while (nodePtr2 != NULL) {
                                    if (nodePtr2->getNameNum() != datum)
                                        nodePtr1->printJacMNAMat(datum, nodePtr2, lastnode, jacMat);
                                    nodePtr2 = nodePtr2->getNext();
                                }
                            }
                            nodePtr1 = nodePtr1->getNext();
                        }
                    }

                    newtonRaphson(nodeValue, number, jacMat, result, delta);        //��ѭ�����������

                    
                }
                for (int j = 0; j < Capacitor_count; j++) {
                    double temp = nodeValue[Ca[0][j]] - nodeValue[Ca[1][j]];
                    Vs[j].push_back(temp);                         //������ͷ�ļ�U(k)
                }
                
            
        }

        //��ӡ���
        for (int i = 0; i < Capacitor_count; i++) {
            cout << "C" << i + 1 << ":" << endl;
            for (int j = 0; j < Vs[i].size(); j++) {
                cout << Vs[i][j] << " ";
            }
            cout << endl;
        }
       

    }



    return 0;
}


// ����Ŀ�꺯������ F(x)
Eigen::VectorXd FX(double result[], int number) {
    Eigen::VectorXd F(number);
    for (int i = 0; i < number; i++) {
        F(i) = result[i + 1];
    }
    return F;
}

// ����Ŀ�꺯�����ſɱȾ��� J(x)��������
Eigen::MatrixXd JA(double jacMat[][30], int number) {
    Eigen::MatrixXd J(number, number);
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < number; j++) {
            J(i, j) = jacMat[i + 1][j + 1];
        }
    }

    return J.inverse();
}

// ţ��-����ɭ��������
void newtonRaphson(double nodeValue[], int number, double jacMat[][30], double result[], double delta[]) {

    Eigen::MatrixXd J = JA(jacMat, number);
    Eigen::VectorXd F = FX(result, number);

    for (int i = 0; i < number; i++) {             //����NR������ϵʽ�������ſɱȾ������ �� F(x)�Ľ��
        for (int j = 0; j < number; j++) {
            delta[i] += J(i, j) * F(j);
        }
    }

    for (int i = 0; i < number; ++i) {
        nodeValue[i + 1] = nodeValue[i + 1] - delta[i];
    }
}

// �жϾ���
bool isAccurate(double delta[], int num, double acc) {
    bool re = true;
    for (int i = 0; i < num; i++) {
        if (delta[i] > acc || -delta[i] > acc) {
            re = false;
        }
    }
    return re;

}

//ͬ�׵���
void homotopy(double nodeValue[], int number, double jacMat[][30], double result[], double delta[], double ratio, const double result_0[]) {

    Eigen::MatrixXd J = JA(jacMat, number);
    Eigen::VectorXd F = FX(result, number);

    for (int i = 0; i < number; i++) {             //����NR������ϵʽ�������ſɱȾ������ �� F(x)�Ľ��
        for (int j = 0; j < number; j++) {
            delta[i] += J(i, j) * (F(j) + ratio * result_0[j + 1]);
        }
    }

    for (int i = 0; i < number; ++i) {
        nodeValue[i + 1] = nodeValue[i + 1] - delta[i];
    }
}

//��˲̬����
void simpleTran(double Uc[], double h, double E, double C, double R, double stop_time, int& size) {
    for (int i = 0; i < (int)(stop_time / h); i++) {
        Uc[i + 1] = h * E / (R * C) + (1 - h / (R * C)) * Uc[i];
        size++;
    }
    /*for (int i = 0; i < (int)(stop_time / h); i++) {
        cout << Uc[i] << endl;
    }*/
}

//%*************************************************************************************************





double stripString(char* stringIn) {
    char buf[BufLength], buf2[BufLength];
    int a, b;
    strcpy(buf, stringIn);
    for (a = 0; buf[a] != '='; a++) {};
    a++;
    for (b = 0; buf[a] != '\0'; b++, a++)
        buf2[b] = buf[a];
    buf2[b] = '\0';
    return atof(buf2);    
};


//Print the linked list of components to check
void printComponents(Component* compPtr) {
    char compTypeName[6];
    cout << endl << "Components: " << endl << endl;
    while (compPtr != NULL) {
        strcpy(compTypeName, strComponentType(compPtr));
        cout << "->" << compTypeName << compPtr->getcompNum();
        compPtr = compPtr->getNext();
    }
    cout << endl;
    return;
}

void printNodes(Node* nodePtr, int compFlag) {

    Connections* conPtr;
    cout << endl << "Nodes: " << endl << endl;
    while (nodePtr != NULL) {
        if (compFlag == 0) { //It is printed just the names of the nodes
            cout << "-> " << nodePtr->getNameNum();
        }
        else if (compFlag == 1) { //It is printed the nodes and the connections
            cout << "-> " << nodePtr->getNameNum() << " {";
            conPtr = nodePtr->getConList();
            while (conPtr->next != NULL) {
                cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << ", ";
                conPtr = conPtr->next;
            }
            cout << strComponentType(conPtr->comp) << conPtr->comp->getcompNum() << '}' << endl;
        }
        else {
            cout << "Invalid value for compFlag. (0) to print just nodes, (1) to print nodes and connections!";
            exit(1);

        }

        nodePtr = nodePtr->getNext();
    }


    return;
}


char* strComponentType(Component* compPtr) {

    char* compTypeName = new char[6];
    switch (compPtr->getType()) {

    case VSource: strcpy(compTypeName, "V"); break;
    case Resistor: strcpy(compTypeName, "R"); break;
    case BJT: strcpy(compTypeName, "T"); break;
    case MOSFET: strcpy(compTypeName, "M"); break;
    case ISource: strcpy(compTypeName, "I"); break;
    case Inductor: strcpy(compTypeName, "ind"); break;
    case Diode: strcpy(compTypeName, "Diode"); break;
    case Capacitor: strcpy(compTypeName, "C"); break;
    }

    return compTypeName;
}



char* ComponentTypeName(Component* compPtr) {

    char* compTypeName = new char[6];
    switch (compPtr->getType()) {

    case VSource: strcpy(compTypeName, "VSource"); break;
    case Resistor: strcpy(compTypeName, "Resistor"); break;
    case BJT: strcpy(compTypeName, "BJT"); break;
    case MOSFET: strcpy(compTypeName, "MOSFET"); break;
    case ISource: strcpy(compTypeName, "ISource"); break;
    case Inductor: strcpy(compTypeName, "Inductor"); break;
    case Diode: strcpy(compTypeName, "Diode"); break;
    case Capacitor: strcpy(compTypeName, "Capacitor"); break;
    }

    return compTypeName;
}



int connectNum(Component* comPtr, Node* nodePtr) {
    if (comPtr->getNodeNum(0) == nodePtr->getNum()) {
        return 0;
    }
    else {
        return 1;
    }
}

#endif