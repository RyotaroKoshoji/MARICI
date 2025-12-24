#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_CRYSTALLINECONSTRAINTMANAGER_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_CRYSTALLINECONSTRAINTMANAGER_H

#include <vector>

#include "ArgumentOutOfRangeException.h"

#include "CrystalStructure.h"

#include "ConstrainerIndices.h"
#include "ConstrainingAtom.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			namespace Internal
			{
				class CrystallineConstraintManager :public ChemToolkit::Crystallography::CrystalStructure<ConstrainingAtom>
				{
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Constructors, destructor, and operators

				protected:
					CrystallineConstraintManager() noexcept;
					explicit CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell&) noexcept;
					CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell&, const std::vector<ConstrainingAtom>&) noexcept;
					CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell&, std::vector<ConstrainingAtom>&&) noexcept;

				public:
					virtual ~CrystallineConstraintManager() = default;

					CrystallineConstraintManager(const CrystallineConstraintManager&) = default;
					CrystallineConstraintManager(CrystallineConstraintManager&&) noexcept = default;
					CrystallineConstraintManager& operator=(const CrystallineConstraintManager&) = default;
					CrystallineConstraintManager& operator=(CrystallineConstraintManager&&) noexcept = default;

				// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Property

					double feasibleErrorRate() const noexcept;
					double exclusiveRadiusRatio() const noexcept;
					double interatomicDistanceTracerCutoffRatio() const noexcept;
					double interatomicDistanceConstrainerCutoffRatio() const noexcept;
					const std::vector<ConstrainerIndices<TranslatedAtomIndex>>& constrainingIndexPairs() const noexcept;

					void setFeasibleErrorRate(const double);
					void setExclusiveRadiusRatio(const double);
					void setInteratomicDistanceTracerCutoffRatio(const double);
					void setInteratomicDistanceConstrainerCutoffRatio(const double);

				// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Methods

					void normalizeAverageFractionalCoordinates() override;
					void normalizeFractionalCoordinates() override;


					void updateTracingIndexPairs();
					void clearInteratomicDistanceConstraints() noexcept;

					bool isIonicAttractive(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isIonicAttractive(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isIonicRepulsive(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isIonicRepulsive(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;

					bool isConstrainableIonicRepulsionDistance(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isConstrainableCovalentExclusionDistance(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isConstrainableIonicExclusionDistance(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;

				// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Protected methods

				protected:
					void updateConstrainingIndexPairs();

					bool isInnateChemicalBondable(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isInnateCovalentBondable(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isInnateIonicBondable(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;

					bool isConstrainableCovalentBondingDistance(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isConstrainableCovalentBondingDistance(const OriginalAtomIndex, const OriginalAtomIndex, const NumericalVector& translationVector) const noexcept;
					bool isConstrainableIonicBondingDistance(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isConstrainableIonicBondingDistance(const OriginalAtomIndex, const OriginalAtomIndex, const NumericalVector& translationVector) const noexcept;
					bool isConstrainableIonicRepulsionDistance(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isConstrainableCovalentExclusionDistance(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isConstrainableIonicExclusionDistance(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;

					bool isFeasibleCovalentBond(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isFeasibleCovalentBond(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isFeasibleCovalentBond(const OriginalAtomIndex, const OriginalAtomIndex, const double errorRate) const noexcept;
					bool isFeasibleCovalentBond(const OriginalAtomIndex, const TranslatedAtomIndex&, const double errorRate) const noexcept;
					bool isFeasibleIonicBond(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isFeasibleIonicBond(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isFeasibleIonicBond(const OriginalAtomIndex, const OriginalAtomIndex, const double errorRate) const noexcept;
					bool isFeasibleIonicBond(const OriginalAtomIndex, const TranslatedAtomIndex&, const double errorRate) const noexcept;
					bool isFeasibleIonicRepulsion(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isFeasibleIonicRepulsion(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isFeasibleCovalentExclusion(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isFeasibleCovalentExclusion(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;
					bool isFeasibleIonicExclusion(const OriginalAtomIndex, const OriginalAtomIndex) const noexcept;
					bool isFeasibleIonicExclusion(const OriginalAtomIndex, const TranslatedAtomIndex&) const noexcept;

					void createCovalentBond(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void createCovalentBond(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;
					void createIonicBond(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void createIonicBond(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;
					void createIonicRepulsion(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void createIonicRepulsion(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;

					void eraseCovalentBond(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void eraseCovalentBond(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;
					void eraseIonicBond(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void eraseIonicBond(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;
					void eraseIonicRepulsion(const OriginalAtomIndex, const OriginalAtomIndex) noexcept;
					void eraseIonicRepulsion(const OriginalAtomIndex, const TranslatedAtomIndex&) noexcept;

					void clearCovalentBonds() noexcept;
					void clearIonicBonds() noexcept;
					void clearIonicRepulsions() noexcept;

				// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Private methods

				private:
					bool isTraceableCovalentExclusionDistance(const OriginalAtomIndex, const OriginalAtomIndex, const LatticePoint& latticePoint) const noexcept;
					bool isTraceableIonicExclusionDistance(const OriginalAtomIndex, const OriginalAtomIndex, const LatticePoint& latticePoint) const noexcept;
					bool isTraceableIonicRepulsionDistance(const OriginalAtomIndex, const OriginalAtomIndex, const LatticePoint& latticePoint) const noexcept;

					std::vector<LatticePoint> enumerateNeighborLatticePoints(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const NumericalMatrix& inverseBasisVectors, const double neighborZoneRadius) const;

				// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

				private:
					double _feasibleErrorRate;
					double _exclusiveRadiusRatio;
					double _interatomicDistanceTracerCutoffRatio;
					double _interatomicDistanceConstrainerCutoffRatio;

					std::vector<ConstrainerIndices<TranslatedAtomIndex>> _constrainingIndexPairs;
					std::vector<ConstrainerIndices<TranslatedAtomIndex>> _tracingIndexPairs;
				};
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
// Property

inline double MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::feasibleErrorRate() const noexcept
{
	return _feasibleErrorRate;
}

inline double MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::exclusiveRadiusRatio() const noexcept
{
	return _exclusiveRadiusRatio;
}

inline double MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::interatomicDistanceTracerCutoffRatio() const noexcept
{
	return _interatomicDistanceTracerCutoffRatio;
}

inline double MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::interatomicDistanceConstrainerCutoffRatio() const noexcept
{
	return _interatomicDistanceConstrainerCutoffRatio;
}

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainerIndices<MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::TranslatedAtomIndex>>& MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::constrainingIndexPairs() const noexcept
{
	return _constrainingIndexPairs;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::setFeasibleErrorRate(const double val)
{
	if (val < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setFeasibleErrorRate", "Argument value is less than zero." };

	if (0.0 < val)
		_feasibleErrorRate = val;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::setExclusiveRadiusRatio(const double val)
{
	if (val < 1.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setExclusiveRadiusRatio", "Argument value is less than one." };

	if (0.0 < val)
		_exclusiveRadiusRatio = val;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::setInteratomicDistanceTracerCutoffRatio(const double val)
{
	if (val < 1.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInteratomicDistanceTracerCutoffRatio", "Argument value is less than one." };

	if (1.0 < val)
		_interatomicDistanceTracerCutoffRatio = val;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::setInteratomicDistanceConstrainerCutoffRatio(const double val)
{
	if (val < 1.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInteratomicDistanceConstrainerCutoffRatio", "Argument value is less than one." };

	if (1.0 < val)
		_interatomicDistanceConstrainerCutoffRatio = val;
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

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::normalizeAverageFractionalCoordinates()
{
	clearInteratomicDistanceConstraints();
	CrystalStructure::normalizeAverageFractionalCoordinates();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::normalizeFractionalCoordinates()
{
	clearInteratomicDistanceConstraints();
	CrystalStructure::normalizeFractionalCoordinates();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::clearInteratomicDistanceConstraints() noexcept
{
	_constrainingIndexPairs.clear();
	_tracingIndexPairs.clear();

	clearCovalentBonds();
	clearIonicBonds();
	clearIonicRepulsions();
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isIonicAttractive(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];


	if (originalAtom.ionicAtomicNumber().formalCharge().isAnion() && translatedAtom.ionicAtomicNumber().formalCharge().isCation())
		return true;
	else if (originalAtom.ionicAtomicNumber().formalCharge().isCation() && translatedAtom.ionicAtomicNumber().formalCharge().isAnion())
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isIonicAttractive(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	return isIonicAttractive(originalAtomIndex, translatedAtomIndex.originalIndex());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isIonicRepulsive(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];


	if (originalAtom.ionicAtomicNumber().formalCharge().isAnion() && translatedAtom.ionicAtomicNumber().formalCharge().isAnion())
		return true;
	else if (originalAtom.ionicAtomicNumber().formalCharge().isCation() && translatedAtom.ionicAtomicNumber().formalCharge().isCation())
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isIonicRepulsive(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	return isIonicRepulsive(originalAtomIndex, translatedAtomIndex.originalIndex());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicRepulsionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableCovalentExclusionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicExclusionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Protected methods

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isInnateChemicalBondable(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	if (atoms()[originalAtomIndex].coordinationConstraints().isInfeasibleAtomicNumber(atoms()[translatedOriginalAtomIndex].ionicAtomicNumber().atomicNumber()))
		return false;

	else
	{
		if (atoms()[translatedOriginalAtomIndex].coordinationConstraints().isInfeasibleAtomicNumber(atoms()[originalAtomIndex].ionicAtomicNumber().atomicNumber()))
			return false;
		else
			return true;
	}
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isInnateCovalentBondable(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	if (0 < atoms()[originalAtomIndex].coordinationConstraints().maxConstrainedCovalentCoordinationNumber())
	{
		if (0 < atoms()[translatedOriginalAtomIndex].coordinationConstraints().maxConstrainedCovalentCoordinationNumber())
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isInnateIonicBondable(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	if (0 < atoms()[originalAtomIndex].coordinationConstraints().maxConstrainedIonicCoordinationNumber())
	{
		if (0 < atoms()[translatedOriginalAtomIndex].coordinationConstraints().maxConstrainedIonicCoordinationNumber())
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableCovalentBondingDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableCovalentBondingDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const NumericalVector& translationVector) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + translationVector - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicBondingDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicBondingDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const NumericalVector& translationVector) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + translationVector - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicRepulsionDistance(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableCovalentExclusionDistance(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isConstrainableIonicExclusionDistance(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double maximumDistance = _interatomicDistanceConstrainerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum());
	double maxmumDistance = (1.0 + _feasibleErrorRate) * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum());
	double maxmumDistance = (1.0 + _feasibleErrorRate) * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const double errorRate) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - errorRate) * (originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum());
	double maxmumDistance = (1.0 + errorRate) * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex, const double errorRate) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - errorRate) * (originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum());
	double maxmumDistance = (1.0 + errorRate) * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum());
	double maxmumDistance = (1.0 + _feasibleErrorRate) * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum());
	double maxmumDistance = (1.0 + _feasibleErrorRate) * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const double errorRate) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - errorRate) * (originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum());
	double maxmumDistance = (1.0 + errorRate) * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex, const double errorRate) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - errorRate) * (originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum());
	double maxmumDistance = (1.0 + errorRate) * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
	{
		if (distanceSquare < (maxmumDistance * maxmumDistance))
			return true;
		else
			return false;
	}

	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - _feasibleErrorRate) * (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentExclusion(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - _feasibleErrorRate) * _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleCovalentExclusion(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - _feasibleErrorRate) * _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicExclusion(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double minmumDistance = (1.0 - _feasibleErrorRate) * _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isFeasibleIonicExclusion(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedAtomIndex.originalIndex()];

	double minmumDistance = (1.0 - _feasibleErrorRate) * _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(translatedAtomIndex.latticePoint()) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if ((minmumDistance * minmumDistance) < distanceSquare)
		return true;
	else
		return false;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createCovalentBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].createCovalentBondWith(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].createCovalentBondWith(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createCovalentBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedAtomIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].createCovalentBondWith(translatedAtomIndex);
	atoms()[translatedAtomIndex.originalIndex()].createCovalentBondWith(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createIonicBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].createIonicBondWith(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].createIonicBondWith(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createIonicBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedAtomIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].createIonicBondWith(translatedAtomIndex);
	atoms()[translatedAtomIndex.originalIndex()].createIonicBondWith(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].createIonicRepulsionWith(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].createIonicRepulsionWith(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::createIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedAtomIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedAtomIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].createIonicRepulsionWith(translatedAtomIndex);
	atoms()[translatedAtomIndex.originalIndex()].createIonicRepulsionWith(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseCovalentBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].eraseCovalentBond(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].eraseCovalentBond(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseCovalentBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].eraseCovalentBond(translatedIndex);
	atoms()[translatedIndex.originalIndex()].eraseCovalentBond(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseIonicBond(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].eraseIonicBond(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].eraseIonicBond(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseIonicBond(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].eraseIonicBond(translatedIndex);
	atoms()[translatedIndex.originalIndex()].eraseIonicBond(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex) noexcept
{
	atoms()[originalAtomIndex].eraseIonicRepulsion(translatedOriginalAtomIndex);
	atoms()[translatedOriginalAtomIndex].eraseIonicRepulsion(originalAtomIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::eraseIonicRepulsion(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& translatedIndex) noexcept
{
	TranslatedAtomIndex reverseIndex{ originalAtomIndex, translatedIndex.latticePoint() };
	reverseIndex.reverseLatticePoint();

	atoms()[originalAtomIndex].eraseIonicRepulsion(translatedIndex);
	atoms()[translatedIndex.originalIndex()].eraseIonicRepulsion(reverseIndex);
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::clearCovalentBonds() noexcept
{
	for (auto& atom : atoms())
		atom.clearCovalentBonds();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::clearIonicBonds() noexcept
{
	for (auto& atom : atoms())
		atom.clearIonicBonds();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::clearIonicRepulsions() noexcept
{
	for (auto& atom : atoms())
		atom.clearIonicRepulsions();
}

// Protected methods
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

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isTraceableCovalentExclusionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const LatticePoint& latticePoint) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceTracerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());
	double distanceSquare = (translatedAtom.cartesianCoordinate() + toTranslationVector(latticePoint) - originalAtom.cartesianCoordinate()).normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isTraceableIonicExclusionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const LatticePoint& latticePoint) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceTracerCutoffRatio * _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());
	double distanceSquare = (translatedAtom.cartesianCoordinate() + toTranslationVector(latticePoint) - originalAtom.cartesianCoordinate()).normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::Internal::CrystallineConstraintManager::isTraceableIonicRepulsionDistance(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex translatedOriginalAtomIndex, const LatticePoint& latticePoint) const noexcept
{
	const ConstrainingAtom& originalAtom = atoms()[originalAtomIndex];
	const ConstrainingAtom& translatedAtom = atoms()[translatedOriginalAtomIndex];

	double maximumDistance = _interatomicDistanceTracerCutoffRatio * (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	NumericalVector displacement = translatedAtom.cartesianCoordinate() + toTranslationVector(latticePoint) - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();


	if (distanceSquare < (maximumDistance * maximumDistance))
		return true;
	else
		return false;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_CRYSTALLINECONSTRAINTMANAGER_H
