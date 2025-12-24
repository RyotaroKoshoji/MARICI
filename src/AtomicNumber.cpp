#include "AtomicNumber.h"

#include <limits>

#include "ArgumentOutOfRangeException.h"

using namespace ChemToolkit::Generic;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

AtomicNumber::AtomicNumber() noexcept
	: _number{ 0 }
{
}

AtomicNumber::AtomicNumber(const size_type val)
	: _number{ val }
{
	if (!(isValid(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Atomic number is out of definitions." };
}

AtomicNumber::AtomicNumber(const std::size_t val)
	: _number{ static_cast<size_type>(val) }
{
	if (std::numeric_limits<size_type>::max() < val)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is too large." };

	if (!(isValid(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Atomic number is out of definitions." };
}

AtomicNumber::AtomicNumber(const int val)
	: _number{ static_cast<size_type>(val) }
{
	if (std::numeric_limits<size_type>::max() < val)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is too large." };

	if (val < 0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument value is less than zero." };

	if (!(isValid(_number)))
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Atomic number is out of definitions." };
}

AtomicNumber::AtomicNumber(const std::string elementSymbol)
	: _number{ ElementSymbol{ elementSymbol }.toAtomicNumber() }
{
}

AtomicNumber::AtomicNumber(const ElementSymbol elementSymbol)
	: _number{ elementSymbol.toAtomicNumber() }
{
}

// Constructors
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Methods

ElementSymbol AtomicNumber::toElementSymbol() const
{
	switch (_number)
	{
	case 0:
		return ElementSymbol{ std::string{ "E" } };
	case 1:
		return ElementSymbol{ std::string{ "H" } };
	case 2:
		return ElementSymbol{ std::string{ "He" } };
	case 3:
		return ElementSymbol{ std::string{ "Li" } };
	case 4:
		return ElementSymbol{ std::string{ "Be" } };
	case 5:
		return ElementSymbol{ std::string{ "B" } };
	case 6:
		return ElementSymbol{ std::string{ "C" } };
	case 7:
		return ElementSymbol{ std::string{ "N" } };
	case 8:
		return ElementSymbol{ std::string{ "O" } };
	case 9:
		return ElementSymbol{ std::string{ "F" } };
	case 10:
		return ElementSymbol{ std::string{ "Ne" } };
	case 11:
		return ElementSymbol{ std::string{ "Na" } };
	case 12:
		return ElementSymbol{ std::string{ "Mg" } };
	case 13:
		return ElementSymbol{ std::string{ "Al" } };
	case 14:
		return ElementSymbol{ std::string{ "Si" } };
	case 15:
		return ElementSymbol{ std::string{ "P" } };
	case 16:
		return ElementSymbol{ std::string{ "S" } };
	case 17:
		return ElementSymbol{ std::string{ "Cl" } };
	case 18:
		return ElementSymbol{ std::string{ "Ar" } };
	case 19:
		return ElementSymbol{ std::string{ "K" } };
	case 20:
		return ElementSymbol{ std::string{ "Ca" } };
	case 21:
		return ElementSymbol{ std::string{ "Sc" } };
	case 22:
		return ElementSymbol{ std::string{ "Ti" } };
	case 23:
		return ElementSymbol{ std::string{ "V" } };
	case 24:
		return ElementSymbol{ std::string{ "Cr" } };
	case 25:
		return ElementSymbol{ std::string{ "Mn" } };
	case 26:
		return ElementSymbol{ std::string{ "Fe" } };
	case 27:
		return ElementSymbol{ std::string{ "Co" } };
	case 28:
		return ElementSymbol{ std::string{ "Ni" } };
	case 29:
		return ElementSymbol{ std::string{ "Cu" } };
	case 30:
		return ElementSymbol{ std::string{ "Zn" } };
	case 31:
		return ElementSymbol{ std::string{ "Ga" } };
	case 32:
		return ElementSymbol{ std::string{ "Ge" } };
	case 33:
		return ElementSymbol{ std::string{ "As" } };
	case 34:
		return ElementSymbol{ std::string{ "Se" } };
	case 35:
		return ElementSymbol{ std::string{ "Br" } };
	case 36:
		return ElementSymbol{ std::string{ "Kr" } };
	case 37:
		return ElementSymbol{ std::string{ "Rb" } };
	case 38:
		return ElementSymbol{ std::string{ "Sr" } };
	case 39:
		return ElementSymbol{ std::string{ "Y" } };
	case 40:
		return ElementSymbol{ std::string{ "Zr" } };
	case 41:
		return ElementSymbol{ std::string{ "Nb" } };
	case 42:
		return ElementSymbol{ std::string{ "Mo" } };
	case 43:
		return ElementSymbol{ std::string{ "Tc" } };
	case 44:
		return ElementSymbol{ std::string{ "Ru" } };
	case 45:
		return ElementSymbol{ std::string{ "Rh" } };
	case 46:
		return ElementSymbol{ std::string{ "Pd" } };
	case 47:
		return ElementSymbol{ std::string{ "Ag" } };
	case 48:
		return ElementSymbol{ std::string{ "Cd" } };
	case 49:
		return ElementSymbol{ std::string{ "In" } };
	case 50:
		return ElementSymbol{ std::string{ "Sn" } };
	case 51:
		return ElementSymbol{ std::string{ "Sb" } };
	case 52:
		return ElementSymbol{ std::string{ "Te" } };
	case 53:
		return ElementSymbol{ std::string{ "I" } };
	case 54:
		return ElementSymbol{ std::string{ "Xe" } };
	case 55:
		return ElementSymbol{ std::string{ "Cs" } };
	case 56:
		return ElementSymbol{ std::string{ "Ba" } };
	case 57:
		return ElementSymbol{ std::string{ "La" } };
	case 58:
		return ElementSymbol{ std::string{ "Ce" } };
	case 59:
		return ElementSymbol{ std::string{ "Pr" } };
	case 60:
		return ElementSymbol{ std::string{ "Nd" } };
	case 61:
		return ElementSymbol{ std::string{ "Pm" } };
	case 62:
		return ElementSymbol{ std::string{ "Sm" } };
	case 63:
		return ElementSymbol{ std::string{ "Eu" } };
	case 64:
		return ElementSymbol{ std::string{ "Gd" } };
	case 65:
		return ElementSymbol{ std::string{ "Tb" } };
	case 66:
		return ElementSymbol{ std::string{ "Dy" } };
	case 67:
		return ElementSymbol{ std::string{ "Ho" } };
	case 68:
		return ElementSymbol{ std::string{ "Er" } };
	case 69:
		return ElementSymbol{ std::string{ "Tm" } };
	case 70:
		return ElementSymbol{ std::string{ "Yb" } };
	case 71:
		return ElementSymbol{ std::string{ "Lu" } };
	case 72:
		return ElementSymbol{ std::string{ "Hf" } };
	case 73:
		return ElementSymbol{ std::string{ "Ta" } };
	case 74:
		return ElementSymbol{ std::string{ "W" } };
	case 75:
		return ElementSymbol{ std::string{ "Re" } };
	case 76:
		return ElementSymbol{ std::string{ "Os" } };
	case 77:
		return ElementSymbol{ std::string{ "Ir" } };
	case 78:
		return ElementSymbol{ std::string{ "Pt" } };
	case 79:
		return ElementSymbol{ std::string{ "Au" } };
	case 80:
		return ElementSymbol{ std::string{ "Hg" } };
	case 81:
		return ElementSymbol{ std::string{ "Tl" } };
	case 82:
		return ElementSymbol{ std::string{ "Pb" } };
	case 83:
		return ElementSymbol{ std::string{ "Bi" } };
	case 84:
		return ElementSymbol{ std::string{ "Po" } };
	case 85:
		return ElementSymbol{ std::string{ "At" } };
	case 86:
		return ElementSymbol{ std::string{ "Rn" } };
	case 87:
		return ElementSymbol{ std::string{ "Fr" } };
	case 88:
		return ElementSymbol{ std::string{ "Ra" } };
	case 89:
		return ElementSymbol{ std::string{ "Ac" } };
	case 90:
		return ElementSymbol{ std::string{ "Th" } };
	case 91:
		return ElementSymbol{ std::string{ "Pa" } };
	case 92:
		return ElementSymbol{ std::string{ "U" } };
	default:
		return ElementSymbol{};
	}
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
