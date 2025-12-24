#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_CRYSTALOPTIMIZER_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_CRYSTALOPTIMIZER_H

#include "NumericalVector.h"
#include "NumericalMatrix.h"

#include "ConstrainerIndices.h"

#include "CrystalDesignRecorder.h"

#include "GeometricalConstraintParameters.h"
#include "StructuralOptimizationParameters.h"

#include "ObjectiveCrystalStructure.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Optimization
		{
			class CrystalOptimizer
			{
				using size_type = std::size_t;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
				using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;

				using SphericalAtom = MathematicalCrystalChemistry::CrystalModel::Components::SphericalAtom;
				using ObjectiveCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure;

				using LatticePoint = ChemToolkit::Crystallography::TranslatedAtomIndex::LatticePoint;
				using OriginalConstrainerIndices = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainerIndices<ChemToolkit::Crystallography::OriginalAtomIndex>;
				using TranslatedConstrainerIndices = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainerIndices<ChemToolkit::Crystallography::TranslatedAtomIndex>;

				using CrystalDesignRecorder = MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				CrystalOptimizer() noexcept;
				CrystalOptimizer(const StructuralOptimizationParameters&, const GeometricalConstraintParameters&);
				virtual ~CrystalOptimizer() = default;

				CrystalOptimizer(const CrystalOptimizer&) = default;
				CrystalOptimizer(CrystalOptimizer&&) noexcept = default;
				CrystalOptimizer& operator=(const CrystalOptimizer&) = default;
				CrystalOptimizer& operator=(CrystalOptimizer&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const StructuralOptimizationParameters& structuralOptimizationParameters() const noexcept;

				void setParameters(const StructuralOptimizationParameters&, const GeometricalConstraintParameters&);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void execute(ObjectiveCrystalStructure&) const;
				void execute(ObjectiveCrystalStructure&, CrystalDesignRecorder&) const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void applyForces(ObjectiveCrystalStructure&) const noexcept;
				void applyPressure(const ObjectiveCrystalStructure&) const noexcept;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private utility

				NumericalVector getTranslationVector(const LatticePoint&, const ObjectiveCrystalStructure&) const noexcept;

				void applyCovalentBondingForce(const OriginalConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyCovalentBondingForce(const TranslatedConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyCovalentExcludingForce(const OriginalConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyCovalentExcludingForce(const TranslatedConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;

				void applyIonicBondingForce(const OriginalConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyIonicBondingForce(const TranslatedConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyIonicExcludingForce(const OriginalConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyIonicExcludingForce(const TranslatedConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyIonicRepulsingForce(const OriginalConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;
				void applyIonicRepulsingForce(const TranslatedConstrainerIndices&, ObjectiveCrystalStructure&) const noexcept;


				void applyForce(SphericalAtom&, SphericalAtom&, const double minimumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, NumericalVector&& translationVector, const double minimumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, const double minimumDistance, const double maximumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, NumericalVector&& translationVector, const double minimumDistance, const double maximumDistance) const noexcept;

				void addUnitCellTransformation(const size_type columnIndex, const NumericalVector& displacement) const noexcept;

			// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				StructuralOptimizationParameters _structuralOptimizationParameters;
				double _exclusiveRadiusRatio;

				mutable NumericalMatrix m_inverseBasisVectors;
				mutable NumericalMatrix m_unitCellTransformation;
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

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::structuralOptimizationParameters() const noexcept
{
	return _structuralOptimizationParameters;
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::setParameters(const StructuralOptimizationParameters& structural, const GeometricalConstraintParameters& geometrical)
{
	_structuralOptimizationParameters = structural;
	_exclusiveRadiusRatio = geometrical.minimumExclusionDistanceRatio();
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
// Private methods

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyPressure(const ObjectiveCrystalStructure& structure) const noexcept
{
	const NumericalMatrix& basisVectors = structure.unitCell().basisVectors();


	m_unitCellTransformation(0, 0) += _structuralOptimizationParameters.pressure() * (basisVectors(2, 1) * basisVectors(1, 2) - basisVectors(1, 1) * basisVectors(2, 2));
	m_unitCellTransformation(1, 0) += _structuralOptimizationParameters.pressure() * (basisVectors(0, 1) * basisVectors(2, 2) - basisVectors(2, 1) * basisVectors(0, 2));
	m_unitCellTransformation(2, 0) += _structuralOptimizationParameters.pressure() * (basisVectors(1, 1) * basisVectors(0, 2) - basisVectors(0, 1) * basisVectors(1, 2));

	m_unitCellTransformation(0, 1) += _structuralOptimizationParameters.pressure() * (basisVectors(2, 2) * basisVectors(1, 0) - basisVectors(1, 2) * basisVectors(2, 0));
	m_unitCellTransformation(1, 1) += _structuralOptimizationParameters.pressure() * (basisVectors(0, 2) * basisVectors(2, 0) - basisVectors(2, 2) * basisVectors(0, 0));
	m_unitCellTransformation(2, 1) += _structuralOptimizationParameters.pressure() * (basisVectors(1, 2) * basisVectors(0, 0) - basisVectors(0, 2) * basisVectors(1, 0));

	m_unitCellTransformation(0, 2) += _structuralOptimizationParameters.pressure() * (basisVectors(2, 0) * basisVectors(1, 1) - basisVectors(1, 0) * basisVectors(2, 1));
	m_unitCellTransformation(1, 2) += _structuralOptimizationParameters.pressure() * (basisVectors(0, 0) * basisVectors(2, 1) - basisVectors(2, 0) * basisVectors(0, 1));
	m_unitCellTransformation(2, 2) += _structuralOptimizationParameters.pressure() * (basisVectors(1, 0) * basisVectors(0, 1) - basisVectors(0, 0) * basisVectors(1, 1));
}

// Private methods
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
// Private utility

inline MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::NumericalVector MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::getTranslationVector(const LatticePoint& latticePoint, const ObjectiveCrystalStructure& structure) const noexcept
{
	NumericalMatrix::const_iterator iter = structure.unitCell().basisVectors().begin();
	double latticePointA = static_cast<double>(latticePoint[0]);
	double latticePointB = static_cast<double>(latticePoint[1]);
	double latticePointC = static_cast<double>(latticePoint[2]);

	NumericalVector translationVector;
	{
		translationVector[0] += (*iter) * latticePointA;
		++iter;
		translationVector[1] += (*iter) * latticePointA;
		++iter;
		translationVector[2] += (*iter) * latticePointA;
		++iter;

		translationVector[0] += (*iter) * latticePointB;
		++iter;
		translationVector[1] += (*iter) * latticePointB;
		++iter;
		translationVector[2] += (*iter) * latticePointB;
		++iter;

		translationVector[0] += (*iter) * latticePointC;
		++iter;
		translationVector[1] += (*iter) * latticePointC;
		++iter;
		translationVector[2] += (*iter) * latticePointC;
	}

	return translationVector;
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyCovalentBondingForce(const OriginalConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum();
	const double maximumDistance = originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum();

	applyForce(originalAtom, translatedAtom, minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyCovalentBondingForce(const TranslatedConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex().originalIndex()];

	const double minimumDistance = originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum();
	const double maximumDistance = originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum();

	applyForce(originalAtom, translatedAtom, getTranslationVector(indices.translatedAtomIndex().latticePoint(), structure), minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyCovalentExcludingForce(const OriginalConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyCovalentExcludingForce(const TranslatedConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex().originalIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	applyForce(originalAtom, translatedAtom, getTranslationVector(indices.translatedAtomIndex().latticePoint(), structure), minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicBondingForce(const OriginalConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum();
	const double maximumDistance = originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum();

	applyForce(originalAtom, translatedAtom, minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicBondingForce(const TranslatedConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex().originalIndex()];

	const double minimumDistance = originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum();
	const double maximumDistance = originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum();

	applyForce(originalAtom, translatedAtom, getTranslationVector(indices.translatedAtomIndex().latticePoint(), structure), minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicExcludingForce(const OriginalConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicExcludingForce(const TranslatedConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex().originalIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	applyForce(originalAtom, translatedAtom, getTranslationVector(indices.translatedAtomIndex().latticePoint(), structure), minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicRepulsingForce(const OriginalConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyIonicRepulsingForce(const TranslatedConstrainerIndices& indices, ObjectiveCrystalStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex().originalIndex()];

	const double minimumDistance = originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum();

	applyForce(originalAtom, translatedAtom, getTranslationVector(indices.translatedAtomIndex().latticePoint(), structure), minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, const double minimumDistance) const noexcept
{
	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;


		addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
		addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
		addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, NumericalVector&& displacement, const double minimumDistance) const noexcept
{
	displacement += (translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate());
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;


		addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
		addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
		addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, const double minimumDistance, const double maximumDistance) const noexcept
{
	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;


		addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
		addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
		addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
	}

	else
	{
		if ((maximumDistance * maximumDistance) < distanceSquare)
		{
			NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

			displacement /= std::sqrt(distanceSquare);
			displacement *= _structuralOptimizationParameters.attractiveForceConstant();
			originalAtom.appliedForce() += displacement;
			translatedAtom.appliedForce() -= displacement;


			addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
			addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
			addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
		}
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, NumericalVector&& displacement, const double minimumDistance, const double maximumDistance) const noexcept
{
	displacement += (translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate());
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;


		addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
		addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
		addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
	}

	else
	{
		if ((maximumDistance * maximumDistance) < distanceSquare)
		{
			NumericalVector fractionalDisplacement = m_inverseBasisVectors * displacement;

			displacement /= std::sqrt(distanceSquare);
			displacement *= _structuralOptimizationParameters.attractiveForceConstant();
			originalAtom.appliedForce() += displacement;
			translatedAtom.appliedForce() -= displacement;


			addUnitCellTransformation(0, (fractionalDisplacement[0] * displacement));
			addUnitCellTransformation(1, (fractionalDisplacement[1] * displacement));
			addUnitCellTransformation(2, (fractionalDisplacement[2] * displacement));
		}
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer::addUnitCellTransformation(const size_type columnIndex, const NumericalVector& displacement) const noexcept
{
	auto basisVectorsTransformationIter = m_unitCellTransformation.begin();
	basisVectorsTransformationIter += (3 * columnIndex);
	{
		auto displacementIter = displacement.begin();

		(*basisVectorsTransformationIter) -= (*displacementIter);
		++basisVectorsTransformationIter;
		++displacementIter;
		(*basisVectorsTransformationIter) -= (*displacementIter);
		++basisVectorsTransformationIter;
		++displacementIter;
		(*basisVectorsTransformationIter) -= (*displacementIter);
	}
}

// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_CRYSTALOPTIMIZER_H
