#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTSDICTIONARY_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTSDICTIONARY_H

#include <vector>
#include <unordered_map>
#include <functional>

#include "StreamReader.h"

#include "AtomicNumber.h"
#include "ChemicalComposition.h"

#include "IonicAtomicNumber.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			class CoordinationConstraintsDictionary
			{
				using size_type = unsigned short;

				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;


				struct KeyHasher
				{
					std::size_t operator()(const IonicAtomicNumber&) const noexcept;
				};

				struct KeyEqual
				{
					bool operator()(const IonicAtomicNumber&, const IonicAtomicNumber&) const noexcept;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			private:
				CoordinationConstraintsDictionary() noexcept = delete;

			public:
				virtual ~CoordinationConstraintsDictionary() = default;

				CoordinationConstraintsDictionary(const CoordinationConstraintsDictionary&) = default;
				CoordinationConstraintsDictionary(CoordinationConstraintsDictionary&&) noexcept = default;
				CoordinationConstraintsDictionary& operator=(const CoordinationConstraintsDictionary&) = default;
				CoordinationConstraintsDictionary& operator=(CoordinationConstraintsDictionary&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				static void initialize() noexcept;
				static void initialize(const System::IO::StreamReader&);

				static bool hasFeasibleCoordnationCompositions(const IonicAtomicNumber&) noexcept;
				static bool hasFeasibleCovalentCoordinationNumbers(const IonicAtomicNumber&) noexcept;
				static bool hasFeasibleIonicCoordinationNumbers(const IonicAtomicNumber&) noexcept;
				static bool hasLowerBoundCoordnationCompositions(const IonicAtomicNumber&) noexcept;

				static std::vector<ChemicalComposition> getFeasibleCoordnationCompositions(const IonicAtomicNumber&) noexcept;
				static std::vector<size_type> getFeasibleCovalentCoordinationNumbers(const IonicAtomicNumber&) noexcept;
				static std::vector<size_type> getFeasibleIonicCoordinationNumbers(const IonicAtomicNumber&) noexcept;
				static std::vector<ChemicalComposition> getLowerBoundCoordnationCompositions(const IonicAtomicNumber&) noexcept;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				static size_type toCoordinationNumber(const std::string&);
				static ChemicalComposition toChemicalComposition(const std::string&);

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				static std::unordered_map<IonicAtomicNumber, std::vector<ChemicalComposition>, KeyHasher, KeyEqual> s_feasibleCompositionDictionary;

				static std::unordered_map<IonicAtomicNumber, std::vector<size_type>, KeyHasher, KeyEqual> s_feasibleCovalentCoordinationNumbers;
				static std::unordered_map<IonicAtomicNumber, std::vector<size_type>, KeyHasher, KeyEqual> s_feasibleIonicCoordinationNumbers;
				static std::unordered_map<IonicAtomicNumber, std::vector<ChemicalComposition>, KeyHasher, KeyEqual> s_lowerBoundCompositionDictionary;
			};



			inline std::size_t CoordinationConstraintsDictionary::KeyHasher::operator()(const IonicAtomicNumber& ian) const noexcept
			{
				return ian.getHashCode();
			}

			inline bool CoordinationConstraintsDictionary::KeyEqual::operator()(const IonicAtomicNumber& iana, const IonicAtomicNumber& ianb) const noexcept
			{
				return (iana == ianb);
			}
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
// Methods

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::hasFeasibleCoordnationCompositions(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	if (s_feasibleCompositionDictionary.find(ionicAtomicNumber) == s_feasibleCompositionDictionary.end())
		return false;
	else
		return true;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::hasFeasibleCovalentCoordinationNumbers(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	if (s_feasibleCovalentCoordinationNumbers.find(ionicAtomicNumber) == s_feasibleCovalentCoordinationNumbers.end())
		return false;
	else
		return true;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::hasFeasibleIonicCoordinationNumbers(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	if (s_feasibleIonicCoordinationNumbers.find(ionicAtomicNumber) == s_feasibleIonicCoordinationNumbers.end())
		return false;
	else
		return true;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::hasLowerBoundCoordnationCompositions(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	if (s_lowerBoundCompositionDictionary.find(ionicAtomicNumber) == s_lowerBoundCompositionDictionary.end())
		return false;
	else
		return true;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::ChemicalComposition> MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::getFeasibleCoordnationCompositions(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_feasibleCompositionDictionary.find(ionicAtomicNumber);

	if (iter == s_feasibleCompositionDictionary.end())
		return std::vector<ChemicalComposition>{};
	else
		return iter->second;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::size_type> MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::getFeasibleCovalentCoordinationNumbers(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_feasibleCovalentCoordinationNumbers.find(ionicAtomicNumber);

	if (iter == s_feasibleCovalentCoordinationNumbers.end())
		return std::vector<size_type>{};
	else
		return iter->second;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::size_type> MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::getFeasibleIonicCoordinationNumbers(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_feasibleIonicCoordinationNumbers.find(ionicAtomicNumber);

	if (iter == s_feasibleIonicCoordinationNumbers.end())
		return std::vector<size_type>{};
	else
		return iter->second;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::ChemicalComposition> MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraintsDictionary::getLowerBoundCoordnationCompositions(const IonicAtomicNumber& ionicAtomicNumber) noexcept
{
	auto iter = s_lowerBoundCompositionDictionary.find(ionicAtomicNumber);

	if (iter == s_lowerBoundCompositionDictionary.end())
		return std::vector<ChemicalComposition>{};
	else
		return iter->second;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_COORDINATIONCONSTRAINTSDICTIONARY_H
