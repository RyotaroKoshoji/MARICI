#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTS_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTS_H

#include <limits>
#include <set>

#include "AtomicNumber.h"
#include "ChemicalComposition.h"

#include "IonicAtomicNumber.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			class CoordinationConstraints
			{
				using size_type = unsigned short;
				using charge_type = short;

				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				CoordinationConstraints() noexcept;
				explicit CoordinationConstraints(const IonicAtomicNumber&);
				CoordinationConstraints(const AtomicNumber, const charge_type formalCharge);

				virtual ~CoordinationConstraints() = default;

				CoordinationConstraints(const CoordinationConstraints&) = default;
				CoordinationConstraints(CoordinationConstraints&&) noexcept = default;
				CoordinationConstraints& operator=(const CoordinationConstraints&) = default;
				CoordinationConstraints& operator=(CoordinationConstraints&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				size_type maxCoordinationNumber() const noexcept;

				size_type maxConstrainedCovalentCoordinationNumber() const noexcept;
				size_type maxConstrainedIonicCoordinationNumber() const noexcept;

				const std::set<ChemicalComposition>& feasibleCompositions() const noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				bool hasFeasibleCompositions() const noexcept;
				bool hasLowerBoundCompositions() const noexcept;
				bool hasFeasibleCovalentCoordinationNumbers() const noexcept;
				bool hasFeasibleIonicCoordinationNumbers() const noexcept;

				bool isInfeasibleAtomicNumber(const AtomicNumber) const noexcept;
				bool isInfeasibleComposition(const ChemicalComposition&) const noexcept;
				bool isInfeasibleCovalentCoordinationNumber(const size_type) const noexcept;
				bool isInfeasibleIonicCoordinationNumber(const size_type) const noexcept;
				bool satisfyLowerBound(const ChemicalComposition&) const;

				ChemicalComposition getClosestFeasibleComposition(const ChemicalComposition&) const;
				ChemicalComposition getClosestLowerBoundComposition(const ChemicalComposition& covalentComposition, const ChemicalComposition& ionicComposition) const;
				size_type getClosestFeasibleCovalentCoordinationNumber(const size_type) const;
				size_type getClosestFeasibleIonicCoordinationNumber(const size_type) const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				std::set<ChemicalComposition> _feasibleCompositionDictionary;

				std::set<ChemicalComposition> _lowerBoundCompositionDictionary;
				std::set<size_type> _feasibleCovalentCoordinationNumbers;
				std::set<size_type> _feasibleIonicCoordinationNumbers;

				std::set<AtomicNumber> _componentAtomicNumbers;
				size_type _maxCoordinationNumber;
				size_type _maxConstrainedCovalentCoordinationNumber;
				size_type _maxConstrainedIonicCoordinationNumber;
			};
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

inline MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::size_type MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::maxCoordinationNumber() const noexcept
{
	return _maxCoordinationNumber;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::size_type MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::maxConstrainedCovalentCoordinationNumber() const noexcept
{
	return _maxConstrainedCovalentCoordinationNumber;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::size_type MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::maxConstrainedIonicCoordinationNumber() const noexcept
{
	return _maxConstrainedIonicCoordinationNumber;
}

inline const std::set<MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::ChemicalComposition>& MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::feasibleCompositions() const noexcept
{
	return _feasibleCompositionDictionary;
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

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::hasFeasibleCompositions() const noexcept
{
	return !(_feasibleCompositionDictionary.empty());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::hasLowerBoundCompositions() const noexcept
{
	return !(_lowerBoundCompositionDictionary.empty());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::hasFeasibleCovalentCoordinationNumbers() const noexcept
{
	return !(_feasibleCovalentCoordinationNumbers.empty());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::hasFeasibleIonicCoordinationNumbers() const noexcept
{
	return !(_feasibleIonicCoordinationNumbers.empty());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::isInfeasibleAtomicNumber(const AtomicNumber an) const noexcept
{
	if (_feasibleCompositionDictionary.empty())
		return false;

	else
	{
		if (_componentAtomicNumbers.find(an) == _componentAtomicNumbers.end())
			return true;
		else
			return false;
	}
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::isInfeasibleComposition(const ChemicalComposition& cc) const noexcept
{
	if (_feasibleCompositionDictionary.empty())
		return false;

	else
	{
		if (_feasibleCompositionDictionary.find(cc) == _feasibleCompositionDictionary.end())
			return true;
		else
			return false;
	}
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::isInfeasibleCovalentCoordinationNumber(const size_type num) const noexcept
{
	if (_feasibleCovalentCoordinationNumbers.empty())
		return false;

	else
	{
		if (_feasibleCovalentCoordinationNumbers.find(num) == _feasibleCovalentCoordinationNumbers.end())
			return true;
		else
			return false;
	}
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::isInfeasibleIonicCoordinationNumber(const size_type num) const noexcept
{
	if (_feasibleIonicCoordinationNumbers.empty())
		return false;
	else
	{
		if (_feasibleIonicCoordinationNumbers.find(num) == _feasibleIonicCoordinationNumbers.end())
			return true;
		else
			return false;
	}
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints::satisfyLowerBound(const ChemicalComposition& cc) const
{
	for (const auto& lowerBound : _lowerBoundCompositionDictionary)
	{
		bool satisfy = true;
		{
			for (const auto& numAndCount : lowerBound)
			{
				if (cc.count(numAndCount.first) < numAndCount.second)
				{
					satisfy = false;
					break;
				}
			}
		}

		if (satisfy)
			return true;
	}

	return false;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTS_H
