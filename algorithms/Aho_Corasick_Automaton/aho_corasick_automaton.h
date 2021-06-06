#pragma once
#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <deque>
#include <map>
using namespace std;
vector<char> TWIG = { 'A','T','G','C' };
vector<string> short_reads_container;
vector<int> randomSeed = { 7,13,17,23 };

int FIND = 0;
class Component {
public:
	int length = 0;
	char twig;
	bool visited = false;
	string shortread = "";
	vector<Component*> arc;
	Component* fail_arc;
	Component() : twig('\0') { arc = {}; fail_arc = {}; }
	Component(char a) : twig(a) { arc = {}; fail_arc = {}; }
	Component* find(char a) {
		for (auto i : arc) {
			if (i->twig == a) return i;
		}
		return nullptr;
	}
}; 
vector<Component*> finalStates;
class automata {
public:

	Component* root;
	automata() {
		root = new Component('$');
	}
	void make_automaton(string short_read) {
		Component* temp = root;
		int read_length = short_read.size();
		//cout << temp->twig << endl;
		for (int i = 0; i < read_length; i++) {
			temp->length = i;
			Component* ptr = temp->find(short_read[i]);
			if (ptr == nullptr) {
				Component* newnode = new Component(short_read[i]);

				// change
				// final state info update
				// shortread : shortread string info, length : read_length (shortread size)
				if (i == read_length - 1) {
					newnode->shortread = short_read;
					newnode->length = read_length;
					//cout << finalStates.size() << endl;
					finalStates.push_back(newnode);
				}

				temp->arc.push_back(newnode);
				temp = newnode;
				//cout << temp->twig << endl;
			}
			else if (ptr != nullptr) {
				//cout << "������" << endl;
				temp = ptr;
			}
		}

	}
	void make_failure_route() {
		Component* temp = root;
		deque<Component*> deck = { temp };
		while (!deck.empty()) {

			Component* ptr = deck.front(); deck.pop_front(); //�θ������

			for (auto i : ptr->arc) { //�ڽĳ�� ��ȸ
				//#1 ���ʳ��κ��� ���丶Ÿ�� �ν� ���� �߽߰� ��Ʈ�� ���ư���.
				if (ptr == temp) {
					ptr->fail_arc = temp;
					for (auto i : ptr->arc) {
						i->fail_arc = temp;
					}
				}
				else {
					//#2 �θ����� ���и�ũ�� �̵��ϰ�, �ڽ� �� ���̺ΰ� ���� ������ ���� ã�´�.
					// �� ã�� ��� �̵��� ��忡���� ���и�ũ�� �̵��Ѵ�.
					Component* fail_moved = ptr->fail_arc;
					Component* needtodecide;
					do {
						needtodecide = fail_moved->find(i->twig);
						fail_moved = fail_moved->fail_arc;
					} while (fail_moved != temp && needtodecide == nullptr);

					if (needtodecide != nullptr) {
						fail_moved = needtodecide;
					}
					i->fail_arc = fail_moved;
					//cout << i->twig<<"#" << fail_moved->twig << "#"<< endl;
				}
				deck.push_back(i);
			}
		}
	}


	void print_automata(Component* temp) {

		for (int i = 0; i < temp->arc.size(); i++) {
			cout << temp->arc[i]->twig << " ";
		} cout << endl;
		for (int i = 0; i < temp->arc.size(); i++) {
			if (temp->arc[i] != nullptr)print_automata(temp->arc[i]);
		}


	}
	void print_node(Component* temp) {
		cout << "index : " << temp->length << endl;
		cout << "string : " << temp->shortread << "ShortRead Length : " << temp->shortread.length() << endl;
		cout << "Twig :" << temp->twig << endl;
		cout << "����� ��� :";
		for (auto i : temp->arc) {
			cout << i->twig << " ";
		}cout << endl;
		cout << "���� ��ũ :";
		cout << temp->fail_arc->twig << endl;


	}

};