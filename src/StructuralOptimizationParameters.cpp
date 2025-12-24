#include "StructuralOptimizationParameters.h"

#include <limits>

#include "InvalidFileException.h"

#include "LengthCasting.h"

using namespace MathematicalCrystalChemistry::Design::Optimization;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

double StructuralOptimizationParameters::s_defaultPressure{ 1.0 };
double StructuralOptimizationParameters::s_defaultAttractiveForceConstant{ 30.0 };
double StructuralOptimizationParameters::s_defaultRepulsiveForceConstant{ -100.0 };
double StructuralOptimizationParameters::s_defaultMaxUnitCellDisplacementFactor{ 0.02 };

StructuralOptimizationParameters::size_type StructuralOptimizationParameters::s_defaultGlobalMaxStructuralOptimizing{ 25 };
double StructuralOptimizationParameters::s_defaultGlobalInitialMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.5) };
double StructuralOptimizationParameters::s_defaultGlobalFinalMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.5) };
double StructuralOptimizationParameters::s_defaultGlobalFeasibleGeometricalConstraintErrorRate{ 1.0 };

StructuralOptimizationParameters::size_type StructuralOptimizationParameters::s_defaultLocalMaxStructuralOptimizing{ 2000 };
double StructuralOptimizationParameters::s_defaultLocalInitialMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.3) };
double StructuralOptimizationParameters::s_defaultLocalFinalMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.05) };
double StructuralOptimizationParameters::s_defaultLocalFeasibleGeometricalConstraintErrorRate{ 0.1 };

StructuralOptimizationParameters::size_type StructuralOptimizationParameters::s_defaultPreciseMaxStructuralOptimizing{ 4000 };
double StructuralOptimizationParameters::s_defaultPreciseInitialMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.1) };
double StructuralOptimizationParameters::s_defaultPreciseFinalMaxAtomicDisplacement{ MathToolkit::UnitConversion::LengthCasting::cast<MathToolkit::UnitConversion::LengthCasting::Unit::Angstrom, MathToolkit::UnitConversion::LengthCasting::Unit::AtomicUnit>(0.005) };
double StructuralOptimizationParameters::s_defaultPreciseFeasibleGeometricalConstraintErrorRate{ 0.05 };


StructuralOptimizationParameters::StructuralOptimizationParameters() noexcept
	: _pressure{ s_defaultPressure }
	, _attractiveForceConstant{ s_defaultAttractiveForceConstant }
	, _repulsiveForceConstant{ s_defaultRepulsiveForceConstant }
	, _maxStructuralOptimizing{ 0 }
	, _initialMaxAtomicDisplacement{ 0.0 }
	, _initialMaxUnitCellDisplacement{ 0.0 }
	, _displacementDecreasingFactor{ 0.0 }
	, _feasibleGeometricalConstraintErrorRate{ 0.0 }
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
// Methods

void StructuralOptimizationParameters::initialize() noexcept
{
	_pressure = s_defaultPressure;
	_attractiveForceConstant = s_defaultAttractiveForceConstant;
	_repulsiveForceConstant = s_defaultRepulsiveForceConstant;

	_maxStructuralOptimizing = 0;
	_initialMaxAtomicDisplacement = 0.0;
	_initialMaxUnitCellDisplacement = 0.0;
	_displacementDecreasingFactor = 0.0;
	_feasibleGeometricalConstraintErrorRate = 0.0;
}

void StructuralOptimizationParameters::initialize(const OptimizationType optimizationType)
{
	_pressure = s_defaultPressure;
	_attractiveForceConstant = s_defaultAttractiveForceConstant;
	_repulsiveForceConstant = s_defaultRepulsiveForceConstant;


	if (optimizationType == OptimizationType::global)
	{
		_maxStructuralOptimizing = s_defaultGlobalMaxStructuralOptimizing;
		_initialMaxAtomicDisplacement = s_defaultGlobalInitialMaxAtomicDisplacement;
		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * s_defaultMaxUnitCellDisplacementFactor;
		_displacementDecreasingFactor = std::pow((s_defaultGlobalFinalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / _maxStructuralOptimizing));
		_feasibleGeometricalConstraintErrorRate = s_defaultGlobalFeasibleGeometricalConstraintErrorRate;
	}

	else if (optimizationType == OptimizationType::local)
	{
		_maxStructuralOptimizing = s_defaultLocalMaxStructuralOptimizing;
		_initialMaxAtomicDisplacement = s_defaultLocalInitialMaxAtomicDisplacement;
		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * s_defaultMaxUnitCellDisplacementFactor;
		_displacementDecreasingFactor = std::pow((s_defaultLocalFinalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / _maxStructuralOptimizing));
		_feasibleGeometricalConstraintErrorRate = s_defaultLocalFeasibleGeometricalConstraintErrorRate;
	}

	else
	{
		_maxStructuralOptimizing = s_defaultPreciseMaxStructuralOptimizing;
		_initialMaxAtomicDisplacement = s_defaultPreciseInitialMaxAtomicDisplacement;
		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * s_defaultMaxUnitCellDisplacementFactor;
		_displacementDecreasingFactor = std::pow((s_defaultPreciseFinalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / _maxStructuralOptimizing));
		_feasibleGeometricalConstraintErrorRate = s_defaultPreciseFeasibleGeometricalConstraintErrorRate;
	}
}

void StructuralOptimizationParameters::initialize(const System::IO::StreamReader& streamReader, const OptimizationType optimizationType)
{
	double scalingMultiplier = 0.0;
	{
		using namespace MathToolkit::UnitConversion::LengthCasting;
		scalingMultiplier = getScalingMultiplier<Unit::Angstrom, Unit::AtomicUnit>();
	}


	initialize();
	double maxUnitCellDisplacementFactor = 0.0;
	{
		if (!(streamReader.readParameter("Pressure", _pressure)))
			_pressure = s_defaultPressure;

		if (!(streamReader.readParameter("Attractive.Force.Constants", _attractiveForceConstant)))
			_attractiveForceConstant = s_defaultAttractiveForceConstant;

		if (!(streamReader.readParameter("Repulsive.Force.Constants", _repulsiveForceConstant)))
			_repulsiveForceConstant = s_defaultRepulsiveForceConstant;

		if (!(streamReader.readParameter("Maximum.Unit.Cell.Displacement.Factor", maxUnitCellDisplacementFactor)))
			maxUnitCellDisplacementFactor = s_defaultMaxUnitCellDisplacementFactor;
	}


	if (optimizationType == OptimizationType::global)
	{
		if (!(streamReader.readParameter("Number.of.Iterative.Balance.Steps", _maxStructuralOptimizing)))
			_maxStructuralOptimizing = s_defaultGlobalMaxStructuralOptimizing;

		{
			if (streamReader.readParameter("Initial.Maximum.Atomic.Displacement", _initialMaxAtomicDisplacement))
				_initialMaxAtomicDisplacement *= scalingMultiplier;
			else
				_initialMaxAtomicDisplacement = s_defaultGlobalInitialMaxAtomicDisplacement;
		}

		if (!(streamReader.readParameter("Feasible.Geometrical.Constraint.Error.Rate", _feasibleGeometricalConstraintErrorRate)))
			_feasibleGeometricalConstraintErrorRate = s_defaultGlobalFeasibleGeometricalConstraintErrorRate;


		double finalMaxAtomicDisplacement = 0.0;
		{
			if (streamReader.readParameter("Final.Maximum.Atomic.Displacement", finalMaxAtomicDisplacement))
				finalMaxAtomicDisplacement *= scalingMultiplier;
			else
				finalMaxAtomicDisplacement = s_defaultGlobalFinalMaxAtomicDisplacement;
		}

		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * maxUnitCellDisplacementFactor;


		if ((0 < _maxStructuralOptimizing) && (0.0 < _initialMaxAtomicDisplacement))
			_displacementDecreasingFactor = std::pow((finalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / static_cast<double>(_maxStructuralOptimizing)));
	}


	else if (optimizationType == OptimizationType::local)
	{
		if (!(streamReader.readParameter("Number.of.Iterative.Balance.Steps", _maxStructuralOptimizing)))
			_maxStructuralOptimizing = s_defaultLocalMaxStructuralOptimizing;

		{
			if (streamReader.readParameter("Initial.Maximum.Atomic.Displacement", _initialMaxAtomicDisplacement))
				_initialMaxAtomicDisplacement *= scalingMultiplier;
			else
				_initialMaxAtomicDisplacement = s_defaultLocalInitialMaxAtomicDisplacement;
		}

		if (!(streamReader.readParameter("Feasible.Geometrical.Constraint.Error.Rate", _feasibleGeometricalConstraintErrorRate)))
			_feasibleGeometricalConstraintErrorRate = s_defaultLocalFeasibleGeometricalConstraintErrorRate;


		double finalMaxAtomicDisplacement = 0.0;
		{
			if (streamReader.readParameter("Final.Maximum.Atomic.Displacement", finalMaxAtomicDisplacement))
				finalMaxAtomicDisplacement *= scalingMultiplier;
			else
				finalMaxAtomicDisplacement = s_defaultLocalFinalMaxAtomicDisplacement;
		}

		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * maxUnitCellDisplacementFactor;


		if ((0 < _maxStructuralOptimizing) && (0.0 < _initialMaxAtomicDisplacement))
			_displacementDecreasingFactor = std::pow((finalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / static_cast<double>(_maxStructuralOptimizing)));
	}


	else
	{
		if (!(streamReader.readParameter("Number.of.Iterative.Balance.Steps", _maxStructuralOptimizing)))
			_maxStructuralOptimizing = s_defaultPreciseMaxStructuralOptimizing;

		{
			if (streamReader.readParameter("Initial.Maximum.Atomic.Displacement", _initialMaxAtomicDisplacement))
				_initialMaxAtomicDisplacement *= scalingMultiplier;
			else
				_initialMaxAtomicDisplacement = s_defaultPreciseInitialMaxAtomicDisplacement;
		}

		if (!(streamReader.readParameter("Feasible.Geometrical.Constraint.Error.Rate", _feasibleGeometricalConstraintErrorRate)))
			_feasibleGeometricalConstraintErrorRate = s_defaultPreciseFeasibleGeometricalConstraintErrorRate;


		double finalMaxAtomicDisplacement = 0.0;
		{
			if (streamReader.readParameter("Final.Maximum.Atomic.Displacement", finalMaxAtomicDisplacement))
				finalMaxAtomicDisplacement *= scalingMultiplier;
			else
				finalMaxAtomicDisplacement = s_defaultPreciseFinalMaxAtomicDisplacement;
		}

		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * maxUnitCellDisplacementFactor;


		if ((0 < _maxStructuralOptimizing) && (0.0 < _initialMaxAtomicDisplacement))
			_displacementDecreasingFactor = std::pow((finalMaxAtomicDisplacement / _initialMaxAtomicDisplacement), (1.0 / static_cast<double>(_maxStructuralOptimizing)));
	}


	validateInitializedValues();
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

void StructuralOptimizationParameters::validateInitializedValues() const
{
	if (_pressure < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Pressure\" is less than zero." };

	if (_attractiveForceConstant < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Attractive.Force.Constants\" is less than zero." };

	if (0.0 < _repulsiveForceConstant)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Repulsive.Force.Constants\" is more than zero." };

	if (_initialMaxAtomicDisplacement < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Initial.Maximum.Atomic.Displacement\" is less than zero." };

	if (_displacementDecreasingFactor < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Final.Maximum.Atomic.Displacement\" is less than zero." };

	if (_initialMaxUnitCellDisplacement < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Maximum.Unit.Cell.Displacement.Factor\" is less than zero." };

	if (1.0 < _displacementDecreasingFactor)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Initial.Maximum.Atomic.Displacement\" is less than \"Final.Maximum.Atomic.Displacement\"." };

	if (_feasibleGeometricalConstraintErrorRate < 0.0)
		throw System::IO::InvalidFileException{ typeid(*this), "validateInitializedValues", "\"Feasible.Geometrical.Constraint.Error.Rate\" is zero." };
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
