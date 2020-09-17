//!  WoodyPlant.
/*!
  This class describes all woody processes and parameters.
*/

#ifndef WOODYPLANT_H
#define WOODYPLANT_H

#include "Plant.h"


class WoodyPlant : public Plant
{
    public:
        //! A constructor.
        /*!
          A more elaborate description of the constructor.
        */
        WoodyPlant(int index, const Input& _input, Landscape* _land);
        WoodyPlant(int _myID, string _PFTname, double _xCor, double _yCor, const Input& _input, Landscape* _land);

        //! A destructor.
        /*!
          A more elaborate description of the destructor.
        */
        virtual ~WoodyPlant();

    protected:

    private:
};

#endif // WOODYPLANT_H
