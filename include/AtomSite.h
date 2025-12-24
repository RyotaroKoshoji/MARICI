#ifndef CHEMTOOLKIT_CRYSTALLOGRAPHY_ATOMSITE_H
#define CHEMTOOLKIT_CRYSTALLOGRAPHY_ATOMSITE_H

#include <cmath>

#include "ArgumentOutOfRangeException.h"
#include "MemberAccessException.h"
#include "NumericalVector.h"

#include "AtomicNumber.h"
#include "ElementSymbol.h"
#include "SiteSymmetrySymbol.h"
#include "WyckoffSymbol.h"


namespace ChemToolkit
{
	namespace Crystallography
	{
		class AtomSite
		{
			using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;

			using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
			using SiteSymmetrySymbol = ChemToolkit::Crystallography::Symmetry::SiteSymmetrySymbol;
			using WyckoffSymbol = ChemToolkit::Crystallography::Symmetry::WyckoffSymbol;

		public:
			using size_type = unsigned short;
			using charge_type = short;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			AtomSite() noexcept;
			virtual ~AtomSite() = default;

			AtomSite(const AtomSite&) = default;
			AtomSite(AtomSite&&) noexcept = default;
			AtomSite& operator=(const AtomSite&) = default;
			AtomSite& operator=(AtomSite&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			AtomicNumber atomicNumber() const noexcept;
			charge_type formalCharge() const noexcept;
			NumericalVector fractionalCoordinate() const noexcept;

			SiteSymmetrySymbol siteSymmetrySymbol() const noexcept;
			WyckoffSymbol wyckoffSymbol() const noexcept;

			double siteOccupancy() const noexcept;
			size_type siteMultiplicity() const;


			void setAtomicNumber(const AtomicNumber) noexcept;
			void setFormalCharge(const charge_type) noexcept;
			void setFractionalCoordinate(const NumericalVector) noexcept;

			void setSiteLabel(const std::string&);
			void setSiteSymmetrySymbol(const SiteSymmetrySymbol&) noexcept;
			void setWyckoffSymbol(const WyckoffSymbol) noexcept;

			void setSiteOccupancy(const double);
			void setSiteMultiplicity(const size_type);

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			std::string siteLabel() const;
			std::string toSiteTypeSymbol() const noexcept;

			bool hasSymmetryInformation() const noexcept;
			void clear() noexcept;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			AtomicNumber _atomicNumber;
			charge_type _formalCharge;
			NumericalVector _fractionalCoordinate;

			size_type _siteIndex;
			SiteSymmetrySymbol _siteSymmetrySymbol;
			WyckoffSymbol _wyckoffSymbol;

			double _siteOccupancy;
			size_type _siteMultiplicity;
		};



		inline bool operator<(const AtomSite& asa, const AtomSite& asb)
		{
			return (asa.siteLabel() < asb.siteLabel());
		}

		inline bool operator>(const AtomSite& asa, const AtomSite& asb)
		{
			return (asb < asa);
		}

		inline bool operator==(const AtomSite& asa, const AtomSite& asb)
		{
			return (asa.siteLabel() == asb.siteLabel());
		}

		inline bool operator!=(const AtomSite& asa, const AtomSite& asb)
		{
			return !(asa == asb);
		}
	}
}

// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Property

inline ChemToolkit::Crystallography::AtomSite::AtomicNumber ChemToolkit::Crystallography::AtomSite::atomicNumber() const noexcept
{
	return _atomicNumber;
}

inline ChemToolkit::Crystallography::AtomSite::charge_type ChemToolkit::Crystallography::AtomSite::formalCharge() const noexcept
{
	return _formalCharge;
}

inline ChemToolkit::Crystallography::AtomSite::NumericalVector ChemToolkit::Crystallography::AtomSite::fractionalCoordinate() const noexcept
{
	return _fractionalCoordinate;
}

inline ChemToolkit::Crystallography::AtomSite::SiteSymmetrySymbol ChemToolkit::Crystallography::AtomSite::siteSymmetrySymbol() const noexcept
{
	return _siteSymmetrySymbol;
}

inline ChemToolkit::Crystallography::AtomSite::WyckoffSymbol ChemToolkit::Crystallography::AtomSite::wyckoffSymbol() const noexcept
{
	return _wyckoffSymbol;
}

inline double ChemToolkit::Crystallography::AtomSite::siteOccupancy() const noexcept
{
	return _siteOccupancy;
}

inline ChemToolkit::Crystallography::AtomSite::size_type ChemToolkit::Crystallography::AtomSite::siteMultiplicity() const
{
	if (_siteMultiplicity == 0)
		throw System::ExceptionServices::MemberAccessException{ typeid(*this), "siteMultiplicity", "Site multiplicity is invalid." };
	else
		return _siteMultiplicity;
}

inline void ChemToolkit::Crystallography::AtomSite::setAtomicNumber(const AtomicNumber number) noexcept
{
	_atomicNumber = number;
}

inline void ChemToolkit::Crystallography::AtomSite::setFormalCharge(const charge_type charge) noexcept
{
	_formalCharge = charge;
}

inline void ChemToolkit::Crystallography::AtomSite::setFractionalCoordinate(const NumericalVector coordinate) noexcept
{
	_fractionalCoordinate = coordinate;
}

inline void ChemToolkit::Crystallography::AtomSite::setSiteSymmetrySymbol(const SiteSymmetrySymbol& symbol) noexcept
{
	_siteSymmetrySymbol = symbol;
}

inline void ChemToolkit::Crystallography::AtomSite::setWyckoffSymbol(const WyckoffSymbol symbol) noexcept
{
	_wyckoffSymbol = symbol;
}

inline void ChemToolkit::Crystallography::AtomSite::setSiteOccupancy(const double occupancy)
{
	if (0.0 < occupancy)
	{
		if (1.0 < occupancy)
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setSiteOccupancy", "Site occupancy is more than one." };
		else
			_siteOccupancy = occupancy;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setSiteOccupancy", "Site occupancy is not more than zero." };
}

inline void ChemToolkit::Crystallography::AtomSite::setSiteMultiplicity(const size_type val)
{
	if (val == 0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setSiteMultiplicity", "Site multiplicity is zero." };
	else
		_siteMultiplicity = val;
}

// Property
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

inline std::string ChemToolkit::Crystallography::AtomSite::siteLabel() const
{
	if (_siteIndex == 0)
		throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "siteLabel", "Site index is zero." };

	else
	{
		std::string siteLabelText;
		{
			siteLabelText += _atomicNumber.toElementSymbol();
			siteLabelText += std::to_string(_siteIndex);
		}

		return siteLabelText;
	}
}

inline std::string ChemToolkit::Crystallography::AtomSite::toSiteTypeSymbol() const noexcept
{
	std::string siteTypeSymbol;
	siteTypeSymbol += _atomicNumber.toElementSymbol();
	{
		if (0 < _formalCharge)
		{
			siteTypeSymbol += std::to_string(_formalCharge);
			siteTypeSymbol += "+";
		}

		else
		{
			if (_formalCharge < 0)
			{
				siteTypeSymbol += std::to_string(std::abs(_formalCharge));
				siteTypeSymbol += "-";
			}
		}
	}

	return siteTypeSymbol;
}

inline bool ChemToolkit::Crystallography::AtomSite::hasSymmetryInformation() const noexcept
{
	if (_siteIndex == 0)
		return false;

	if (!(_siteSymmetrySymbol.isValid()))
		return false;

	if (!(_wyckoffSymbol.isValid()))
		return false;

	if (_siteMultiplicity == 0)
		return false;


	return true;
}

inline void ChemToolkit::Crystallography::AtomSite::clear() noexcept
{
	_atomicNumber = 0;
	_formalCharge = 0;
	_fractionalCoordinate = 0.0;

	_siteIndex = 0;
	_siteSymmetrySymbol.set();
	_wyckoffSymbol.set();

	_siteOccupancy = 0.0;
	_siteMultiplicity = 0;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !CHEMTOOLKIT_CRYSTALLOGRAPHY_ATOMSITE_H
