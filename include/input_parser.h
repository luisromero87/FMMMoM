using namespace std;

#ifndef INPUT_PARSER
#define INPUT_PARSER

class InputParser {
public:
	InputParser (int &argc, char **argv);
	string getCmdOption(const string &option) const;
	bool cmdOptionExists(const string &option) const;
private:
	vector <string> tokens;
};

#endif
