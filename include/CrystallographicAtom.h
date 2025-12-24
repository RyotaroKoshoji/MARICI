#ifndef CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICATOM_H
#define CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICATOM_H


#include "ArgumentOutOfRangeException.h"
#include "InvalidOperationException.h"

#include "NumericalVector.h"

#include "IonicAtomicNumber.h"
#include "SiteSymmetrySymbol.h"
#include "WyckoffSymbol.h"


namespace ChemToolkit
{
	namespace Crystallography
	{
		class CrystallographicAtom
		{
		protected:
			using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;

			using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
			using SiteSymmetrySymbol = ChemToolkit::Crystallography::Symmetry::SiteSymmetrySymbol;
			using WyckoffSymbol = ChemToolkit::Crystallography::Symmetry::WyckoffSymbol;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystallographicAtom() noexcept;
			explicit CrystallographicAtom(const IonicAtomicNumber&) noexcept;
			CrystallographicAtom(const IonicAtomicNumber&, const NumericalVector&) noexcept;

			virtual ~CrystallographicAtom() = default;

			CrystallographicAtom(const CrystallographicAtom&) = default;
			CrystallographicAtom(CrystallographicAtom&&) noexcept = default;
			CrystallographicAtom& operator=(const CrystallographicAtom&) = default;
			CrystallographicAtom& operator=(CrystallographicAtom&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			IonicAtomicNumber ionicAtomicNumber() const noexcept;
			const NumericalVector& cartesianCoordinate() const noexcept;
			NumericalVector& cartesianCoordinate() noexcept;

			const std::string& siteLabel() const noexcept;
			const SiteSymmetrySymbol& siteSymmetrySymbol() const noexcept;
			WyckoffSymbol wyckoffSymbol() const noexcept;


			void setSiteLabel(const std::string&) noexcept;
			void setSiteSymmetrySymbol(const SiteSymmetrySymbol&) noexcept;
			void setWyckoffSymbol(const WyckoffSymbol&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			void move(const NumericalVector& displacement) noexcept;

			bool hasSymmetryInformation() const noexcept;
			void clearSymmetryInformation() noexcept;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			IonicAtomicNumber _ionicAtomicNumber;
			NumericalVector _cartesianCoordinate;

			std::string _siteLabel;
			SiteSymmetrySymbol _siteSymmetrySymbol;
			WyckoffSymbol _wyckoffSymbol;
		};



		inline bool operator<(const CrystallographicAtom& caa, const CrystallographicAtom& cab)
		{
			if (!(caa.hasSymmetryInformation()) || !(cab.hasSymmetryInformation()))
				throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Could not compare static atoms since atoms are invalid." };

			else
			{
				if (caa.ionicAtomicNumber() < cab.ionicAtomicNumber())
				{
					if (caa.siteLabel() == cab.siteLabel())
						throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Site labels are same insplite of different atomic numbers." };
					else
						return true;
				}

				else if (cab.ionicAtomicNumber() < caa.ionicAtomicNumber())
				{
					if (caa.siteLabel() == cab.siteLabel())
						throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Site labels are same insplite of different atomic numbers." };
					else
						return false;
				}

				else
				{
					if (caa.siteSymmetrySymbol() < cab.siteSymmetrySymbol())
					{
						if (caa.siteLabel() == cab.siteLabel())
							throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Site labels are the same inspite of different site symmetry." };
						else
							return true;
					}

					else if (cab.siteSymmetrySymbol() < caa.siteSymmetrySymbol())
					{
						if (caa.siteLabel() == cab.siteLabel())
							throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Site labels are the same inspite of different site symmetry." };
						else
							return false;
					}

					else
						return false;
				}
			}
		}

		inline bool operator>(const CrystallographicAtom& caa, const CrystallographicAtom& cab)
		{
			return (cab < caa);
		}

		inline bool operator==(const CrystallographicAtom& caa, const CrystallographicAtom& cab)
		{
			if (!(caa.hasSymmetryInformation()) || !(cab.hasSymmetryInformation()))
				throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator==", "Could not equal static atoms since atoms are invalid." };

			else
			{
				if (caa.ionicAtomicNumber() != cab.ionicAtomicNumber())
				{
					if (caa.siteLabel() == cab.siteLabel())
						throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator==", "Site labels are same insplite of different atomic numbers." };
					else
						return false;
				}

				else
				{
					if (caa.siteLabel() == cab.siteLabel())
					{
						if ((caa.siteSymmetrySymbol() == cab.siteSymmetrySymbol()) && (caa.wyckoffSymbol() == cab.wyckoffSymbol()))
							return true;
						else
							throw System::ExceptionServices::InvalidOperationException{ "ChemToolkit::Crystallography::operator<", "Could not equal static atoms." };
					}

					else
						return false;
				}
			}
		}

		inline bool operator!=(const CrystallographicAtom& caa, const CrystallographicAtom& cab)
		{
			return !(caa == cab);
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
// Constructors

inline ChemToolkit::Crystallography::CrystallographicAtom::CrystallographicAtom() noexcept
	: _ionicAtomicNumber{}
	, _cartesianCoordinate{ 0.0,0.0,0.0 }
	, _siteLabel{}
	, _siteSymmetrySymbol{}
	, _wyckoffSymbol{}
{
}

inline ChemToolkit::Crystallography::CrystallographicAtom::CrystallographicAtom(const IonicAtomicNumber& number) noexcept
	: _ionicAtomicNumber{ number }
	, _cartesianCoordinate{ 0.0,0.0,0.0 }
	, _siteLabel{}
	, _siteSymmetrySymbol{}
	, _wyckoffSymbol{}
{
}

inline ChemToolkit::Crystallography::CrystallographicAtom::CrystallographicAtom(const IonicAtomicNumber& number, const NumericalVector& coordinate) noexcept
	: _ionicAtomicNumber{ number }
	, _cartesianCoordinate{ coordinate }
	, _siteLabel{}
	, _siteSymmetrySymbol{}
	, _wyckoffSymbol{}
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
// Property

inline ChemToolkit::Crystallography::CrystallographicAtom::IonicAtomicNumber ChemToolkit::Crystallography::CrystallographicAtom::ionicAtomicNumber() const noexcept
{
	return _ionicAtomicNumber;
}

inline const ChemToolkit::Crystallography::CrystallographicAtom::NumericalVector& ChemToolkit::Crystallography::CrystallographicAtom::cartesianCoordinate() const noexcept
{
	return _cartesianCoordinate;
}

inline ChemToolkit::Crystallography::CrystallographicAtom::NumericalVector& ChemToolkit::Crystallography::CrystallographicAtom::cartesianCoordinate() noexcept
{
	return _cartesianCoordinate;
}

inline const std::string& ChemToolkit::Crystallography::CrystallographicAtom::siteLabel() const noexcept
{
	return _siteLabel;
}

inline const ChemToolkit::Crystallography::CrystallographicAtom::SiteSymmetrySymbol& ChemToolkit::Crystallography::CrystallographicAtom::siteSymmetrySymbol() const noexcept
{
	return _siteSymmetrySymbol;
}

inline ChemToolkit::Crystallography::CrystallographicAtom::WyckoffSymbol ChemToolkit::Crystallography::CrystallographicAtom::wyckoffSymbol() const noexcept
{
	return _wyckoffSymbol;
}

inline void ChemToolkit::Crystallography::CrystallographicAtom::setSiteLabel(const std::string& siteLabel) noexcept
{
	_siteLabel = siteLabel;
}

inline void ChemToolkit::Crystallography::CrystallographicAtom::setSiteSymmetrySymbol(const SiteSymmetrySymbol& siteSymmetrySymbol) noexcept
{
	_siteSymmetrySymbol = siteSymmetrySymbol;
}

inline void ChemToolkit::Crystallography::CrystallographicAtom::setWyckoffSymbol(const WyckoffSymbol& wyckoffSymbol) noexcept
{
	_wyckoffSymbol = wyckoffSymbol;
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

inline void ChemToolkit::Crystallography::CrystallographicAtom::move(const NumericalVector& displacement) noexcept
{
	_cartesianCoordinate += displacement;

	_siteLabel.clear();
	_siteSymmetrySymbol.set();
	_wyckoffSymbol.set();
}

inline bool ChemToolkit::Crystallography::CrystallographicAtom::hasSymmetryInformation() const noexcept
{
	if (_siteLabel.empty())
		return false;

	if (!(_siteSymmetrySymbol.isValid()))
		return false;

	if (!(_wyckoffSymbol.isValid()))
		return false;


	return true;
}

inline void ChemToolkit::Crystallography::CrystallographicAtom::clearSymmetryInformation() noexcept
{
	_siteLabel.clear();
	_siteSymmetrySymbol.set();
	_wyckoffSymbol.set();
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICATOM_H
