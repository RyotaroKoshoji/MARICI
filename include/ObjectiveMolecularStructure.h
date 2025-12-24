#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVEMOLECULARSTRUCTURE_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVEMOLECULARSTRUCTURE_H

#include <vector>

#include "AtomicNumber.h"
#include "CrystalStructure.h"

#include "SphericalAtom.h"
#include "ConstrainerIndices.h"
#include "CoordinationConstraints.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			class ConstrainingMolecularStructure;


			class ObjectiveMolecularStructure :public ChemToolkit::Crystallography::CrystalStructure<SphericalAtom>
			{
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
				using UnitCell = ChemToolkit::Crystallography::UnitCell;
				using CoordinationConstraints = MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints;

				using OriginalConstrainerIndices = ConstrainerIndices<OriginalAtomIndex>;
				using TranslatedConstrainerIndices = ConstrainerIndices<TranslatedAtomIndex>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				ObjectiveMolecularStructure() noexcept;
				ObjectiveMolecularStructure(const UnitCell&, const std::vector<SphericalAtom>&, const std::vector<IonicAtomicNumber>&, const std::vector<CoordinationConstraints>&) noexcept;
				ObjectiveMolecularStructure(const UnitCell&, std::vector<SphericalAtom>&&, std::vector<IonicAtomicNumber>&&, std::vector<CoordinationConstraints>&&) noexcept;

				explicit ObjectiveMolecularStructure(const ConstrainingMolecularStructure&);

				virtual ~ObjectiveMolecularStructure() = default;

				ObjectiveMolecularStructure(const ObjectiveMolecularStructure&) = default;
				ObjectiveMolecularStructure(ObjectiveMolecularStructure&&) noexcept = default;
				ObjectiveMolecularStructure& operator=(const ObjectiveMolecularStructure&) = default;
				ObjectiveMolecularStructure& operator=(ObjectiveMolecularStructure&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const std::vector<IonicAtomicNumber>& correspondingIonicAtomicNumbers() const noexcept;
				const std::vector<CoordinationConstraints>& correspondingCoordinationConstraints() const noexcept;

				const std::vector<OriginalConstrainerIndices>& covalentBondedIndices() const noexcept;
				const std::vector<OriginalConstrainerIndices>& covalentExcludedIndices() const noexcept;
				const std::vector<OriginalConstrainerIndices>& ionicBondedIndices() const noexcept;
				const std::vector<OriginalConstrainerIndices>& ionicExcludedIndices() const noexcept;
				const std::vector<OriginalConstrainerIndices>& ionicRepulsedIndices() const noexcept;


				void setCovalentBondedIndices(const std::vector<OriginalConstrainerIndices>&) noexcept;
				void setCovalentExcludedIndices(const std::vector<OriginalConstrainerIndices>&) noexcept;
				void setIonicBondedIndices(const std::vector<OriginalConstrainerIndices>&) noexcept;
				void setIonicExcludedIndices(const std::vector<OriginalConstrainerIndices>&) noexcept;
				void setIonicRepulsedIndices(const std::vector<OriginalConstrainerIndices>&) noexcept;

				void setCovalentBondedIndices(std::vector<OriginalConstrainerIndices>&&) noexcept;
				void setCovalentExcludedIndices(std::vector<OriginalConstrainerIndices>&&) noexcept;
				void setIonicBondedIndices(std::vector<OriginalConstrainerIndices>&&) noexcept;
				void setIonicExcludedIndices(std::vector<OriginalConstrainerIndices>&&) noexcept;
				void setIonicRepulsedIndices(std::vector<OriginalConstrainerIndices>&&) noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				bool isValid() const noexcept;
				bool isFeasible(const double feasibleErrorRate, const double exclusionRatio) const;

				void import(const ConstrainingMolecularStructure&);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				bool isFeasibleDistance(const OriginalAtomIndex, const OriginalAtomIndex, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept;
				bool isFeasibleDistance(const OriginalAtomIndex, const OriginalAtomIndex, const double minimumDistance, const double feasibleErrorRate) const noexcept;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				std::vector<IonicAtomicNumber> _correspondingIonicAtomicNumbers;
				std::vector<CoordinationConstraints> _correspondingCoordinationConstraints;

				std::vector<OriginalConstrainerIndices> _covalentBondedIndices;
				std::vector<OriginalConstrainerIndices> _covalentExcludedIndices;
				std::vector<OriginalConstrainerIndices> _ionicBondedIndices;
				std::vector<OriginalConstrainerIndices> _ionicExcludedIndices;
				std::vector<OriginalConstrainerIndices> _ionicRepulsedIndices;
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

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::IonicAtomicNumber>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::correspondingIonicAtomicNumbers() const noexcept
{
	return _correspondingIonicAtomicNumbers;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::CoordinationConstraints>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::correspondingCoordinationConstraints() const noexcept
{
	return _correspondingCoordinationConstraints;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::covalentBondedIndices() const noexcept
{
	return _covalentBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::covalentExcludedIndices() const noexcept
{
	return _covalentExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::ionicBondedIndices() const noexcept
{
	return _ionicBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::ionicExcludedIndices() const noexcept
{
	return _ionicExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::ionicRepulsedIndices() const noexcept
{
	return _ionicRepulsedIndices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setCovalentBondedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_covalentBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setCovalentExcludedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_covalentExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicBondedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicExcludedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicRepulsedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicRepulsedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setCovalentBondedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_covalentBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setCovalentExcludedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_covalentExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicBondedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicExcludedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::setIonicRepulsedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicRepulsedIndices = std::move(indices);
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::isValid() const noexcept
{
	if (atoms().size() != _correspondingIonicAtomicNumbers.size())
		return false;

	return true;
}

// Methods
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
// Private methods

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept
{
	NumericalVector displacement = atoms()[translatedIndex].cartesianCoordinate() - atoms()[originalIndex].cartesianCoordinate();
	double distanceSquare = displacement.normSquare();

	double minimumDistanceSquare = (1 - feasibleErrorRate) * (1 - feasibleErrorRate) * minimumDistance * minimumDistance;
	double maximumDistanceSquare = (1 + feasibleErrorRate) * (1 + feasibleErrorRate) * maximumDistance * maximumDistance;


	if (distanceSquare < minimumDistanceSquare)
		return false;
	else if (maximumDistanceSquare < distanceSquare)
		return false;
	else
		return true;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const double minimumDistance, const double feasibleErrorRate) const noexcept
{
	NumericalVector displacement = atoms()[translatedIndex].cartesianCoordinate() - atoms()[originalIndex].cartesianCoordinate();

	double distanceSquare = displacement.normSquare();
	double minimumDistanceSquare = (1 - feasibleErrorRate) * (1 - feasibleErrorRate) * minimumDistance * minimumDistance;


	if (distanceSquare < minimumDistanceSquare)
		return false;
	else
		return true;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVEMOLECULARSTRUCTURE_H
