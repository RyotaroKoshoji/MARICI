#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_GENERATION_RANDOMSTRUCTUREGENERATIONPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_GENERATION_RANDOMSTRUCTUREGENERATIONPARAMETERS_H

#include <filesystem>

#include "ArgumentOutOfRangeException.h"

#include "StreamReader.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Generation
		{
			class RandomStructureGenerationParameters
			{
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				RandomStructureGenerationParameters() noexcept;
				virtual ~RandomStructureGenerationParameters() = default;

				RandomStructureGenerationParameters(const RandomStructureGenerationParameters&) = default;
				RandomStructureGenerationParameters(RandomStructureGenerationParameters&&) noexcept = default;
				RandomStructureGenerationParameters& operator=(const RandomStructureGenerationParameters&) = default;
				RandomStructureGenerationParameters& operator=(RandomStructureGenerationParameters&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				double initialPackingFraction() const noexcept;
				double maxInitialLatticeLengthRatio() const noexcept;
				double maxInitialLatticeAngle() const noexcept;
				double minInitialLatticeAngle() const noexcept;

				void setInitialPackingFraction(const double);
				void setMaxInitialLatticeLengthRatio(const double);
				void setInitialLatticeAngles(const double minAngle, const double maxAngle);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void initialize() noexcept;
				void initialize(const System::IO::StreamReader&);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void validateParameters() const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				double _initialPackingFraction;
				double _maxInitialLatticeLengthRatio;
				double _maxInitialLatticeAngle;
				double _minInitialLatticeAngle;

				static double s_defaultInitialPackingFraction;
				static double s_defaultMaxInitialLatticeLengthRatio;
				static double s_defaultMaxInitialLatticeAngle;
				static double s_defaultMinInitialLatticeAngle;
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

inline double MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::initialPackingFraction() const noexcept
{
	return _initialPackingFraction;
}

inline double MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::maxInitialLatticeLengthRatio() const noexcept
{
	return _maxInitialLatticeLengthRatio;
}

inline double MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::maxInitialLatticeAngle() const noexcept
{
	return _maxInitialLatticeAngle;
}

inline double MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::minInitialLatticeAngle() const noexcept
{
	return _minInitialLatticeAngle;
}

inline void MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::setInitialPackingFraction(const double val)
{
	if (0.0 < val)
	{
		if (val < 1.0)
			_initialPackingFraction = val;
		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInitialPackingFraction", "Initial.Packing.Fraction is not less than one." };
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInitialPackingFraction", "Initial.Packing.Fraction is not more than zero." };
}

inline void MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerationParameters::setMaxInitialLatticeLengthRatio(const double val)
{
	if (_maxInitialLatticeLengthRatio < 1.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMaxInitialLatticeLengthRatio", "Maximum.Initial.Lattice.Length.Ratio is less than one." };
	else
		_maxInitialLatticeLengthRatio = val;
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_GENERATION_RANDOMSTRUCTUREGENERATIONPARAMETERS_H
