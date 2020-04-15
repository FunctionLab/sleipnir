/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef COLOR_H
#define COLOR_H

#include "colori.h"

namespace Sleipnir {

/*!
 * \brief
 * Simple representation of a color triple in RGB space.
 */
    class CColor : CColorImpl {
    public:
        /*!
         * \brief
         * Black
         */
        static const CColor c_Black;
        /*!
         * \brief
         * Cyan
         */
        static const CColor c_Cyan;
        /*!
         * \brief
         * Green
         */
        static const CColor c_Green;
        /*!
         * \brief
         * Red
         */
        static const CColor c_Red;
        /*!
         * \brief
         * White
         */
        static const CColor c_White;
        /*!
         * \brief
         * Yellow
         */
        static const CColor c_Yellow;

        static const CColor c_Blue;

        static const CColor c_DarkGreen;

        static const CColor c_Orange;

        static CColor Interpolate(float dValue, const CColor &ColorMinimum, const CColor &ColorMedium,
                                  const CColor &ColorMaximum);

        CColor(const unsigned char *abRGB);

        CColor(unsigned char bRed, unsigned char bGreen, unsigned char bBlue);

        CColor operator+(const CColor &Color) const;

        CColor operator*(float dValue) const;

        CColor &operator=(const CColor &Color);

        std::string ToRGB() const;

        bool IsDark() const;
    };

}

#endif // COLOR_H
