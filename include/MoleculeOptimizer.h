#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_MOLECULEOPTIMIZER_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_MOLECULEOPTIMIZER_H

#include "NumericalVector.h"
#include "NumericalMatrix.h"

#include "ConstrainerIndices.h"

#include "CrystalDesignRecorder.h"

#include "GeometricalConstraintParameters.h"
#include "StructuralOptimizationParameters.h"

#include "ObjectiveMolecularStructure.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Optimization
		{
			class MoleculeOptimizer
			{
				using size_type = std::size_t;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;
				using NumericalMatrix = MathToolkit::LinearAlgebra::NumericalMatrix<double, 3, 3>;

				using SphericalAtom = MathematicalCrystalChemistry::CrystalModel::Components::SphericalAtom;
				using ObjectiveMolecularStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure;

				using LatticePoint = ChemToolkit::Crystallography::TranslatedAtomIndex::LatticePoint;
				using OriginalConstrainerIndices = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainerIndices<ChemToolkit::Crystallography::OriginalAtomIndex>;

				using CrystalDesignRecorder = MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				MoleculeOptimizer() noexcept;
				MoleculeOptimizer(const StructuralOptimizationParameters&, const GeometricalConstraintParameters&);
				virtual ~MoleculeOptimizer() = default;

				MoleculeOptimizer(const MoleculeOptimizer&) = default;
				MoleculeOptimizer(MoleculeOptimizer&&) noexcept = default;
				MoleculeOptimizer& operator=(const MoleculeOptimizer&) = default;
				MoleculeOptimizer& operator=(MoleculeOptimizer&&) noexcept = default;

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

				void execute(ObjectiveMolecularStructure&) const;
				void execute(ObjectiveMolecularStructure&, CrystalDesignRecorder&) const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void applyForces(ObjectiveMolecularStructure&) const noexcept;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private utility

				NumericalVector getTranslationVector(const LatticePoint&, const ObjectiveMolecularStructure&) const noexcept;

				void applyCovalentBondingForce(const OriginalConstrainerIndices&, ObjectiveMolecularStructure&) const noexcept;
				void applyCovalentExcludingForce(const OriginalConstrainerIndices&, ObjectiveMolecularStructure&) const noexcept;
				void applyIonicBondingForce(const OriginalConstrainerIndices&, ObjectiveMolecularStructure&) const noexcept;
				void applyIonicExcludingForce(const OriginalConstrainerIndices&, ObjectiveMolecularStructure&) const noexcept;
				void applyIonicRepulsingForce(const OriginalConstrainerIndices&, ObjectiveMolecularStructure&) const noexcept;


				void applyForce(SphericalAtom&, SphericalAtom&, const double minimumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, NumericalVector&& translationVector, const double minimumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, const double minimumDistance, const double maximumDistance) const noexcept;
				void applyForce(SphericalAtom&, SphericalAtom&, NumericalVector&& translationVector, const double minimumDistance, const double maximumDistance) const noexcept;

			// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				StructuralOptimizationParameters _structuralOptimizationParameters;
				double _exclusiveRadiusRatio;

				mutable NumericalMatrix m_inverseBasisVectors;
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

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::structuralOptimizationParameters() const noexcept
{
	return _structuralOptimizationParameters;
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::setParameters(const StructuralOptimizationParameters& structural, const GeometricalConstraintParameters& geometrical)
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
// Private utility

inline MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::NumericalVector MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::getTranslationVector(const LatticePoint& latticePoint, const ObjectiveMolecularStructure& structure) const noexcept
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

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyCovalentBondingForce(const OriginalConstrainerIndices& indices, ObjectiveMolecularStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = originalAtom.covalentRadius().minimum() + translatedAtom.covalentRadius().minimum();
	const double maximumDistance = originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum();

	applyForce(originalAtom, translatedAtom, minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyCovalentExcludingForce(const OriginalConstrainerIndices& indices, ObjectiveMolecularStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.covalentRadius().maximum() + translatedAtom.covalentRadius().maximum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyIonicBondingForce(const OriginalConstrainerIndices& indices, ObjectiveMolecularStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = originalAtom.ionicRadius().minimum() + translatedAtom.ionicRadius().minimum();
	const double maximumDistance = originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum();

	applyForce(originalAtom, translatedAtom, minimumDistance, maximumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyIonicExcludingForce(const OriginalConstrainerIndices& indices, ObjectiveMolecularStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = _exclusiveRadiusRatio * (originalAtom.ionicRadius().maximum() + translatedAtom.ionicRadius().maximum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyIonicRepulsingForce(const OriginalConstrainerIndices& indices, ObjectiveMolecularStructure& structure) const noexcept
{
	SphericalAtom& originalAtom = structure.atoms()[indices.originalAtomIndex()];
	SphericalAtom& translatedAtom = structure.atoms()[indices.translatedAtomIndex()];

	const double minimumDistance = (originalAtom.ionicRepulsionRadius().minimum() + translatedAtom.ionicRepulsionRadius().minimum());

	applyForce(originalAtom, translatedAtom, minimumDistance);
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, const double minimumDistance) const noexcept
{
	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, NumericalVector&& displacement, const double minimumDistance) const noexcept
{
	displacement += (translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate());
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, const double minimumDistance, const double maximumDistance) const noexcept
{
	NumericalVector displacement = translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate();
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;
	}

	else
	{
		if ((maximumDistance * maximumDistance) < distanceSquare)
		{
			displacement /= std::sqrt(distanceSquare);
			displacement *= _structuralOptimizationParameters.attractiveForceConstant();
			originalAtom.appliedForce() += displacement;
			translatedAtom.appliedForce() -= displacement;
		}
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::MoleculeOptimizer::applyForce(SphericalAtom& originalAtom, SphericalAtom& translatedAtom, NumericalVector&& displacement, const double minimumDistance, const double maximumDistance) const noexcept
{
	displacement += (translatedAtom.cartesianCoordinate() - originalAtom.cartesianCoordinate());
	double distanceSquare = displacement.normSquare();

	if (distanceSquare < (minimumDistance * minimumDistance))
	{
		displacement /= std::sqrt(distanceSquare);
		displacement *= _structuralOptimizationParameters.repulsiveForceConstant();
		originalAtom.appliedForce() += displacement;
		translatedAtom.appliedForce() -= displacement;
	}

	else
	{
		if ((maximumDistance * maximumDistance) < distanceSquare)
		{
			displacement /= std::sqrt(distanceSquare);
			displacement *= _structuralOptimizationParameters.attractiveForceConstant();
			originalAtom.appliedForce() += displacement;
			translatedAtom.appliedForce() -= displacement;
		}
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


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_MOLECULEOPTIMIZER_H
