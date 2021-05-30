#include <iostream>
#include <fstream>
#include <string>

using namespace std;
int main(void)
{

	int count = 0;
	int mis = 0;
	ifstream p;
	string name = "reference.txt";
	string ref;
	p.open(name);
	p >> ref;
	for (int i = 0; i < 100000000; i++) {
		if (ref[i] == 'A' || ref[i] == 'T' || ref[i] == 'G' || ref[i] == 'C') {
			count++;
		}
		else { mis++; }
	}

	cout << "count : " << count << "mis : " << mis << endl;
	return 0;
}