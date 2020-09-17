#include <iostream>
#include "Input.h"
#include "Landscape.h"

using namespace std;

/////////////////////////////////////////////////////////
int main()
{
    //cout << "Goodday Australia!" << endl;

    //Read in everything
    Input* input;
    input = new Input();

    //Create landscape
    Landscape* land;
    land = new Landscape(*input);

    //Run simulation
    for(int day = 0; day < input->simDays; day++) {
        cout << "day " << day+1 << " of " << input->simDays << endl;
        land->calculateProcesses(day);
    }

    //Delete objects
    delete land;
    land = NULL;
    delete input;
    input = NULL;

    return 0;
}
