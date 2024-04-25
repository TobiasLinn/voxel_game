#include <iostream>

#include "game.hpp"

using namespace std;
using namespace vxg;

int main(int argc, char **argv) {
    try {
        game g;
        g.run();
    } catch (const string &msg) {
        cout << msg << endl;
    }
    return 0;
}
