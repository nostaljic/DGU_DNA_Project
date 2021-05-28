#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(void) {

	ifstream p;
	ofstream f;
	string file = "chr9.fna";
	string file2 = "reference.txt";

	string reference;
	p.open(file); 	// Reference.txt open
	f.open(file2);
	int count = 0;
	for (int i = 0; i < 1403501; i++) {
		p >> reference;
		cout << "count : " << i << reference << endl;
		for (int j = 0; j < 80; j++) {
			if (reference[j] == 'A' || reference[j] == 'T' || reference[j] == 'G' || reference[j] == 'C') {
				f << reference[j];
				//cout << reference[j];
				count++;
				if (count == 100000000) {
					cout << "finish";
					f.close();
					p.close();
					return 0; }
			}
		}
	}
	return 0;
}