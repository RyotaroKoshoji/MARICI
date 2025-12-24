#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVECRYSTALSTRUCTURE_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVECRYSTALSTRUCTURE_H

#include <vector>

#include "IonicAtomicNumber.h"
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
			class ConstrainingCrystalStructure;


			class ObjectiveCrystalStructure :public ChemToolkit::Crystallography::CrystalStructure<SphericalAtom>
			{
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;

				using UnitCell = ChemToolkit::Crystallography::UnitCell;
				using CoordinationConstraints = MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints;

				using OriginalConstrainerIndices = ConstrainerIndices<OriginalAtomIndex>;
				using TranslatedConstrainerIndices = ConstrainerIndices<TranslatedAtomIndex>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				ObjectiveCrystalStructure() noexcept;
				ObjectiveCrystalStructure(const UnitCell&, const std::vector<SphericalAtom>&, const std::vector<IonicAtomicNumber>&, const std::vector<CoordinationConstraints>&) noexcept;
				ObjectiveCrystalStructure(const UnitCell&, std::vector<SphericalAtom>&&, std::vector<IonicAtomicNumber>&&, std::vector<CoordinationConstraints>&&) noexcept;

				explicit ObjectiveCrystalStructure(const ConstrainingCrystalStructure&);

				virtual ~ObjectiveCrystalStructure() = default;

				ObjectiveCrystalStructure(const ObjectiveCrystalStructure&) = default;
				ObjectiveCrystalStructure(ObjectiveCrystalStructure&&) noexcept = default;
				ObjectiveCrystalStructure& operator=(const ObjectiveCrystalStructure&) = default;
				ObjectiveCrystalStructure& operator=(ObjectiveCrystalStructure&&) noexcept = default;

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

				const std::vector<TranslatedConstrainerIndices>& translatedCovalentBondedIndices() const noexcept;
				const std::vector<TranslatedConstrainerIndices>& translatedCovalentExcludedIndices() const noexcept;
				const std::vector<TranslatedConstrainerIndices>& translatedIonicBondedIndices() const noexcept;
				const std::vector<TranslatedConstrainerIndices>& translatedIonicExcludedIndices() const noexcept;
				const std::vector<TranslatedConstrainerIndices>& translatedIonicRepulsedIndices() const noexcept;

				
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

				void setTranslatedCovalentBondedIndices(const std::vector<TranslatedConstrainerIndices>&) noexcept;
				void setTranslatedCovalentExcludedIndices(const std::vector<TranslatedConstrainerIndices>&) noexcept;
				void setTranslatedIonicBondedIndices(const std::vector<TranslatedConstrainerIndices>&) noexcept;
				void setTranslatedIonicExcludedIndices(const std::vector<TranslatedConstrainerIndices>&) noexcept;
				void setTranslatedIonicRepulsedIndices(const std::vector<TranslatedConstrainerIndices>&) noexcept;

				void setTranslatedCovalentBondedIndices(std::vector<TranslatedConstrainerIndices>&&) noexcept;
				void setTranslatedCovalentExcludedIndices(std::vector<TranslatedConstrainerIndices>&&) noexcept;
				void setTranslatedIonicBondedIndices(std::vector<TranslatedConstrainerIndices>&&) noexcept;
				void setTranslatedIonicExcludedIndices(std::vector<TranslatedConstrainerIndices>&&) noexcept;
				void setTranslatedIonicRepulsedIndices(std::vector<TranslatedConstrainerIndices>&&) noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				bool isValid() const noexcept;
				bool isFeasible(const double feasibleErrorRate, const double exclusionRatio) const;

				void import(const ConstrainingCrystalStructure&);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				bool isFeasibleDistance(const OriginalAtomIndex, const OriginalAtomIndex, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept;
				bool isFeasibleDistance(const OriginalAtomIndex, const TranslatedAtomIndex&, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept;
				bool isFeasibleDistance(const OriginalAtomIndex, const OriginalAtomIndex, const double minimumDistance, const double feasibleErrorRate) const noexcept;
				bool isFeasibleDistance(const OriginalAtomIndex, const TranslatedAtomIndex&, const double minimumDistance, const double feasibleErrorRate) const noexcept;

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

				std::vector<TranslatedConstrainerIndices> _translatedCovalentBondedIndices;
				std::vector<TranslatedConstrainerIndices> _translatedCovalentExcludedIndices;
				std::vector<TranslatedConstrainerIndices> _translatedIonicBondedIndices;
				std::vector<TranslatedConstrainerIndices> _translatedIonicExcludedIndices;
				std::vector<TranslatedConstrainerIndices> _translatedIonicRepulsedIndices;
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

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::IonicAtomicNumber>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::correspondingIonicAtomicNumbers() const noexcept
{
	return _correspondingIonicAtomicNumbers;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::CoordinationConstraints>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::correspondingCoordinationConstraints() const noexcept
{
	return _correspondingCoordinationConstraints;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::covalentBondedIndices() const noexcept
{
	return _covalentBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::covalentExcludedIndices() const noexcept
{
	return _covalentExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::ionicBondedIndices() const noexcept
{
	return _ionicBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::ionicExcludedIndices() const noexcept
{
	return _ionicExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::OriginalConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::ionicRepulsedIndices() const noexcept
{
	return _ionicRepulsedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::TranslatedConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::translatedCovalentBondedIndices() const noexcept
{
	return _translatedCovalentBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::TranslatedConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::translatedCovalentExcludedIndices() const noexcept
{
	return _translatedCovalentExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::TranslatedConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::translatedIonicBondedIndices() const noexcept
{
	return _translatedIonicBondedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::TranslatedConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::translatedIonicExcludedIndices() const noexcept
{
	return _translatedIonicExcludedIndices;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::TranslatedConstrainerIndices>& MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::translatedIonicRepulsedIndices() const noexcept
{
	return _translatedIonicRepulsedIndices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setCovalentBondedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_covalentBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setCovalentExcludedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_covalentExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicBondedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicExcludedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicRepulsedIndices(const std::vector<OriginalConstrainerIndices>& indices) noexcept
{
	_ionicRepulsedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setCovalentBondedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_covalentBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setCovalentExcludedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_covalentExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicBondedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicExcludedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setIonicRepulsedIndices(std::vector<OriginalConstrainerIndices>&& indices) noexcept
{
	_ionicRepulsedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedCovalentBondedIndices(const std::vector<TranslatedConstrainerIndices>& indices) noexcept
{
	_translatedCovalentBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedCovalentExcludedIndices(const std::vector<TranslatedConstrainerIndices>& indices) noexcept
{
	_translatedCovalentExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicBondedIndices(const std::vector<TranslatedConstrainerIndices>& indices) noexcept
{
	_translatedIonicBondedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicExcludedIndices(const std::vector<TranslatedConstrainerIndices>& indices) noexcept
{
	_translatedIonicExcludedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicRepulsedIndices(const std::vector<TranslatedConstrainerIndices>& indices) noexcept
{
	_translatedIonicRepulsedIndices = indices;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedCovalentBondedIndices(std::vector<TranslatedConstrainerIndices>&& indices) noexcept
{
	_translatedCovalentBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedCovalentExcludedIndices(std::vector<TranslatedConstrainerIndices>&& indices) noexcept
{
	_translatedCovalentExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicBondedIndices(std::vector<TranslatedConstrainerIndices>&& indices) noexcept
{
	_translatedIonicBondedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicExcludedIndices(std::vector<TranslatedConstrainerIndices>&& indices) noexcept
{
	_translatedIonicExcludedIndices = std::move(indices);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::setTranslatedIonicRepulsedIndices(std::vector<TranslatedConstrainerIndices>&& indices) noexcept
{
	_translatedIonicRepulsedIndices = std::move(indices);
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::isValid() const noexcept
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const TranslatedAtomIndex& translatedIndex, const double minimumDistance, const double maximumDistance, const double feasibleErrorRate) const noexcept
{
	NumericalVector displacement = atoms()[translatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(translatedIndex.latticePoint()) - atoms()[originalIndex].cartesianCoordinate();
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const double minimumDistance, const double feasibleErrorRate) const noexcept
{
	NumericalVector displacement = atoms()[translatedIndex].cartesianCoordinate() - atoms()[originalIndex].cartesianCoordinate();

	double distanceSquare = displacement.normSquare();
	double minimumDistanceSquare = (1 - feasibleErrorRate) * (1 - feasibleErrorRate) * minimumDistance * minimumDistance;


	if (distanceSquare < minimumDistanceSquare)
		return false;
	else
		return true;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure::isFeasibleDistance(const OriginalAtomIndex originalIndex, const TranslatedAtomIndex& translatedIndex, const double minimumDistance, const double feasibleErrorRate) const noexcept
{
	NumericalVector displacement = atoms()[translatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(translatedIndex.latticePoint()) - atoms()[originalIndex].cartesianCoordinate();

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


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OBJECTIVECRYSTALSTRUCTURE_H
