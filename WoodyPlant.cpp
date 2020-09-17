#include "WoodyPlant.h"
#include "Input.h"

/////////////////////////////////////////////////////////
WoodyPlant::WoodyPlant(int index, const Input& _input, Landscape* _land) : Plant(index, _input, _land)
{
    //ctor
}

/////////////////////////////////////////////////////////
WoodyPlant::WoodyPlant(int _myID, string _PFTname, double _xCor, double _yCor, const Input& _input, Landscape* _land) : Plant(_myID, _PFTname, _xCor, _yCor, _input, _land)
{
    //ctor
}

/////////////////////////////////////////////////////////
WoodyPlant::~WoodyPlant()
{
    //dtor
}
