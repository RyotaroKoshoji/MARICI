#ifndef CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICSTRUCTURE_H
#define CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICSTRUCTURE_H

#include <vector>

#include "ElementSymbol.h"
#include "ChemicalComposition.h"

#include "SpaceGroupNumber.h"
#include "SpaceGroupTypeSearcher.h"

#include "UnitCell.h"
#include "CrystallographicInformation.h"


namespace ChemToolkit
{
	namespace Crystallography
	{
		template <typename A>
		class CrystallographicStructure
		{
			using Atom = A;

		protected:
			using size_type = std::size_t;
			using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
			using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;

			using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
			using ElementSymbol = ChemToolkit::Generic::ElementSymbol;
			using SpaceGroupNumber = ChemToolkit::Crystallography::Symmetry::SpaceGroupNumber;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystallographicStructure() noexcept;
			explicit CrystallographicStructure(const UnitCell&) noexcept;
			CrystallographicStructure(const UnitCell&, const std::vector<Atom>&) noexcept;
			CrystallographicStructure(const UnitCell&, std::vector<Atom>&&) noexcept;

			template <typename AA>
			CrystallographicStructure(const CrystallographicStructure<AA>&);


			virtual ~CrystallographicStructure() = default;

			CrystallographicStructure(const CrystallographicStructure&) = default;
			CrystallographicStructure(CrystallographicStructure&&) noexcept = default;
			CrystallographicStructure& operator=(const CrystallographicStructure&) = default;
			CrystallographicStructure& operator=(CrystallographicStructure&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			const UnitCell& unitCell() const noexcept;
			const std::vector<Atom>& atoms() const noexcept;

			bool hasSymmetryInformation() const noexcept;
			SpaceGroupNumber spaceGroupNumber() const noexcept;


			void setUnitCell(const UnitCell&) noexcept;
			void setAtoms(const std::vector<Atom>&) noexcept;
			void setAtoms(std::vector<Atom>&&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			virtual void initialize() noexcept;


			virtual void applyDelaunayReduction();
			virtual void applyDelaunayReduction(const double spaceGroupPrecision);
			virtual void simplifyStructure();
			virtual void simplifyStructure(const double spaceGroupPrecision);
			virtual void conventionalizeStructure();
			virtual void conventionalizeStructure(const double spaceGroupPrecision);
			virtual void updateSymmetryInformation();
			virtual void updateSymmetryInformation(const double spaceGroupPrecision);
			void clearSymmetryInformation() noexcept;


			virtual void normalizeFractionalCoordinates();
			CrystallographicInformation toCrystallographicInformation() const;
			CrystallographicInformation toCrystallographicInformation(const double spaceGroupPrecision) const;

			ChemToolkit::Generic::ChemicalComposition<AtomicNumber> toChemicalComposition() const;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Utility

			virtual void eraseAtoms(const AtomicNumber&);
			virtual void eraseAtoms(const ElementSymbol&);

		// Utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			UnitCell _unitCell;
			std::vector<Atom> _atoms;

			bool _hasSymmetryInformation;
			SpaceGroupNumber _spaceGroupNumber;
		};
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

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicStructure<A>::CrystallographicStructure() noexcept
	: _unitCell{}
	, _atoms{}
	, _hasSymmetryInformation{ false }
	, _spaceGroupNumber{}
{
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicStructure<A>::CrystallographicStructure(const UnitCell& cell) noexcept
	: _unitCell{ cell }
	, _atoms{}
	, _hasSymmetryInformation{ false }
	, _spaceGroupNumber{}
{
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicStructure<A>::CrystallographicStructure(const UnitCell& cell, const std::vector<Atom>& arrangement) noexcept
	: _unitCell{ cell }
	, _atoms{ arrangement }
	, _hasSymmetryInformation{ false }
	, _spaceGroupNumber{}
{
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicStructure<A>::CrystallographicStructure(const UnitCell& cell, std::vector<Atom>&& arrangement) noexcept
	: _unitCell{ cell }
	, _atoms{ std::move(arrangement) }
	, _hasSymmetryInformation{ false }
	, _spaceGroupNumber{}
{
}

template <typename A>
template <typename B>
inline ChemToolkit::Crystallography::CrystallographicStructure<A>::CrystallographicStructure(const CrystallographicStructure<B>& crystalStructure)
	: _unitCell{ crystalStructure.unitCell() }
	, _atoms{}
	, _hasSymmetryInformation{ crystalStructure._hasSymmetryInformation }
	, _spaceGroupNumber{ crystalStructure._spaceGroupNumber }
{
	if (_hasSymmetryInformation)
	{
		for (const auto& inputAtom : crystalStructure.atoms())
		{
			if (inputAtom.hasSymmetryInformation())
			{
				Atom atom{ inputAtom.ionicAtomicNumber(), inputAtom.cartesianCoordinate() };
				{
					atom.setSiteLabel(inputAtom.siteLabel());
					atom.setSiteSymmetrySymbol(inputAtom.siteSymmetrySymbol());
					atom.setWyckoffSymbol(inputAtom.wyckoffSymbol());
				}

				_atoms.push_back(atom);
			}

			else
				throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "Argument atoms do not have the symmetry information." };
		}
	}

	else
	{
		for (const auto& inputAtom : crystalStructure.atoms())
			_atoms.push_back(Atom{ inputAtom.ionicAtomicNumber(), inputAtom.cartesianCoordinate() });
	}
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
// Access methods

template <typename A>
inline const ChemToolkit::Crystallography::UnitCell& ChemToolkit::Crystallography::CrystallographicStructure<A>::unitCell() const noexcept
{
	return _unitCell;
}

template <typename A>
inline const typename std::vector<A>& ChemToolkit::Crystallography::CrystallographicStructure<A>::atoms() const noexcept
{
	return _atoms;
}

template <typename A>
inline bool ChemToolkit::Crystallography::CrystallographicStructure<A>::hasSymmetryInformation() const noexcept
{
	return _hasSymmetryInformation;
}

template <typename A>
inline typename ChemToolkit::Crystallography::CrystallographicStructure<A>::SpaceGroupNumber ChemToolkit::Crystallography::CrystallographicStructure<A>::spaceGroupNumber() const noexcept
{
	return _spaceGroupNumber;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::setUnitCell(const UnitCell& cell) noexcept
{
	_unitCell = cell;
	
	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::setAtoms(const std::vector<Atom>& arrangement) noexcept
{
	_atoms = arrangement;

	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::setAtoms(std::vector<Atom>&& arrangement) noexcept
{
	_atoms = std::move(arrangement);

	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

// Access methods
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
// Public methods

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::initialize() noexcept
{
	_unitCell.clear();
	_atoms.clear();
	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::applyDelaunayReduction()
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher;
	spaceGroupTypeSearcher.applyDelaunayReduction(_unitCell, _atoms);

	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::applyDelaunayReduction(const double spaceGroupPrecision)
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ spaceGroupPrecision };
	spaceGroupTypeSearcher.applyDelaunayReduction(_unitCell, _atoms);

	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::simplifyStructure()
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher;
	spaceGroupTypeSearcher.toPrimitive(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::simplifyStructure(const double spaceGroupPrecision)
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ spaceGroupPrecision };
	spaceGroupTypeSearcher.toPrimitive(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::conventionalizeStructure()
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher;
	spaceGroupTypeSearcher.conventionalize(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::conventionalizeStructure(const double spaceGroupPrecision)
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ spaceGroupPrecision };
	spaceGroupTypeSearcher.conventionalize(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::updateSymmetryInformation()
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher;
	spaceGroupTypeSearcher.updateSymmetryInformation(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::updateSymmetryInformation(const double spaceGroupPrecision)
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ spaceGroupPrecision };
	spaceGroupTypeSearcher.updateSymmetryInformation(_unitCell, _atoms, _spaceGroupNumber);

	_hasSymmetryInformation = true;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::clearSymmetryInformation() noexcept
{
	for (auto& atom : _atoms)
		atom.clearSymmetryInformation();


	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::normalizeFractionalCoordinates()
{
	for (auto& atom : _atoms)
	{
		NumericalVector fractionalCoordinate = _unitCell.getInverseBasisVectors() * atom.cartesianCoordinate();
		int floorValueA = static_cast<int>(std::floor(fractionalCoordinate[0]));
		int floorValueB = static_cast<int>(std::floor(fractionalCoordinate[1]));
		int floorValueC = static_cast<int>(std::floor(fractionalCoordinate[2]));

		if (floorValueA != 0)
			fractionalCoordinate[0] -= static_cast<double>(floorValueA);
		if (floorValueB != 0)
			fractionalCoordinate[1] -= static_cast<double>(floorValueB);
		if (floorValueC != 0)
			fractionalCoordinate[2] -= static_cast<double>(floorValueC);


		atom.cartesianCoordinate() = (_unitCell.basisVectors() * fractionalCoordinate);
	}


	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicInformation ChemToolkit::Crystallography::CrystallographicStructure<A>::toCrystallographicInformation() const
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher;

	if (_hasSymmetryInformation)
		return spaceGroupTypeSearcher.getCrystallographicInformation(_unitCell, _atoms, _spaceGroupNumber);
	else
		return spaceGroupTypeSearcher.getCrystallographicInformation(_unitCell, _atoms);
}

template <typename A>
inline ChemToolkit::Crystallography::CrystallographicInformation ChemToolkit::Crystallography::CrystallographicStructure<A>::toCrystallographicInformation(const double spaceGroupPrecision) const
{
	ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ spaceGroupPrecision };

	if (_hasSymmetryInformation)
		return spaceGroupTypeSearcher.getCrystallographicInformation(_unitCell, _atoms, _spaceGroupNumber);
	else
		return spaceGroupTypeSearcher.getCrystallographicInformation(_unitCell, _atoms);
}

template <typename A>
inline ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> ChemToolkit::Crystallography::CrystallographicStructure<A>::toChemicalComposition() const
{
	ChemToolkit::Generic::ChemicalComposition<AtomicNumber> chemicalComposition;
	{
		for (const auto& atom : _atoms)
			chemicalComposition.add(atom.ionicAtomicNumber().atomicNumber(), 1);
	}

	return chemicalComposition;
}

// Public methods
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
// Utility

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::eraseAtoms(const AtomicNumber& number)
{
	for (auto iter = _atoms.begin(); iter != _atoms.end();)
	{
		if (iter->ionicAtomicNumber().atomicNumber() == number)
			iter = _atoms.erase(iter);
		else
			++iter;
	}


	_hasSymmetryInformation = false;
	_spaceGroupNumber = 1;
}

template <typename A>
inline void ChemToolkit::Crystallography::CrystallographicStructure<A>::eraseAtoms(const ElementSymbol& symbol)
{
	eraseAtoms(symbol.toAtomicNumber());
}

// Utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !CHEMTOOLKIT_CRYSTALLOGRAPHY_CRYSTALLOGRAPHICSTRUCTURE_H
