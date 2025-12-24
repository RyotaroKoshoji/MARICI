#ifndef MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_PROMISINGCRYSTALEXTRACTIONPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_PROMISINGCRYSTALEXTRACTIONPARAMETERS_H

#include "StreamReader.h"

#include "SpaceGroupNumber.h"

#include <utility>


namespace MathematicalCrystalChemistry
{
	namespace Analysis
	{
		class PromisingCrystalExtractionParameters
		{
			using size_type = unsigned short;
			using SpaceGroupNumber = ChemToolkit::Crystallography::Symmetry::SpaceGroupNumber;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			PromisingCrystalExtractionParameters() noexcept;
			virtual ~PromisingCrystalExtractionParameters() = default;

			PromisingCrystalExtractionParameters(const PromisingCrystalExtractionParameters&) = default;
			PromisingCrystalExtractionParameters(PromisingCrystalExtractionParameters&&) noexcept = default;
			PromisingCrystalExtractionParameters& operator=(const PromisingCrystalExtractionParameters&) = default;
			PromisingCrystalExtractionParameters& operator=(PromisingCrystalExtractionParameters&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			bool needExtraction() const noexcept;

			bool identifySiteLabel() const noexcept;
			bool identifySiteSymmetry() const noexcept;
			bool identifyCoordinations() const noexcept;
			bool identifyCoordinationPolyhedraLinkings() const noexcept;

			double spaceGroupPrecision() const noexcept;
			const std::pair<SpaceGroupNumber, SpaceGroupNumber>& feasibleSpaceGroupRange() const noexcept;

			size_type maximumAtomicEnvironments() const noexcept;
			size_type maximumAnionicEnvironments() const noexcept;
			size_type maximumCationicEnvironments() const noexcept;

			size_type maximumAtomicEnvironmentsPerElement() const noexcept;
			size_type maximumAtomicEnvironmentsPerAnionicElement() const noexcept;
			size_type maximumAtomicEnvironmentsPerCationicElement() const noexcept;


			void setSpaceGroupPrecision(const double);
			void setSiteLabelIdentification(const bool) noexcept;
			void setSiteSymmetryIdentification(const bool) noexcept;
			void setCoordinationIdentification(const bool) noexcept;
			void setCoordinationPolyhedraLinkingsIdentification(const bool) noexcept;

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
			bool toNecessity(const std::string&) const;
			void validate() const;

		// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			bool _needExtraction;

			bool _identifySiteLabel;
			bool _identifySiteSymmetry;
			bool _identifyCoordinations;
			bool _identifyCoordinationPolyhedraLinkings;

			double _spaceGroupPrecision;
			std::pair<SpaceGroupNumber, SpaceGroupNumber> _feasibleSpaceGroupRange;

			size_type _maximumAtomicEnvironments;
			size_type _maximumAnionicEnvironments;
			size_type _maximumCationicEnvironments;

			size_type _maximumAtomicEnvironmentsPerElement;
			size_type _maximumAtomicEnvironmentsPerAnionicElement;
			size_type _maximumAtomicEnvironmentsPerCationicElement;


			static double s_defaultSpaceGroupPrecision;
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
// Property

inline bool MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::needExtraction() const noexcept
{
	return _needExtraction;
}

inline bool MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::identifySiteLabel() const noexcept
{
	return _identifySiteLabel;
}

inline bool MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::identifySiteSymmetry() const noexcept
{
	return _identifySiteSymmetry;
}

inline bool MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::identifyCoordinations() const noexcept
{
	return _identifyCoordinations;
}

inline bool MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::identifyCoordinationPolyhedraLinkings() const noexcept
{
	return _identifyCoordinationPolyhedraLinkings;
}

inline double MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::spaceGroupPrecision() const noexcept
{
	return _spaceGroupPrecision;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumAtomicEnvironments() const noexcept
{
	return _maximumAtomicEnvironments;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumAnionicEnvironments() const noexcept
{
	return _maximumAnionicEnvironments;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumCationicEnvironments() const noexcept
{
	return _maximumCationicEnvironments;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumAtomicEnvironmentsPerElement() const noexcept
{
	return _maximumAtomicEnvironmentsPerElement;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumAtomicEnvironmentsPerAnionicElement() const noexcept
{
	return _maximumAtomicEnvironmentsPerAnionicElement;
}

inline MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::size_type MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::maximumAtomicEnvironmentsPerCationicElement() const noexcept
{
	return _maximumAtomicEnvironmentsPerCationicElement;
}

inline const std::pair<MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::SpaceGroupNumber, MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::SpaceGroupNumber>& MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::feasibleSpaceGroupRange() const noexcept
{
	return _feasibleSpaceGroupRange;
}

inline void MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::setSiteLabelIdentification(const bool necessity) noexcept
{
	_identifySiteLabel = necessity;
}

inline void MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::setSiteSymmetryIdentification(const bool necessity) noexcept
{
	_identifySiteSymmetry = necessity;
}

inline void MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::setCoordinationIdentification(const bool necessity) noexcept
{
	_identifyCoordinations = necessity;
}

inline void MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::setCoordinationPolyhedraLinkingsIdentification(const bool necessity) noexcept
{
	_identifyCoordinationPolyhedraLinkings = necessity;
}

inline void MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters::setSpaceGroupPrecision(const double val)
{
	if (0.0 < val)
		_spaceGroupPrecision = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setSpaceGroupPrecision", "Argument value is not more than zero." };
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_PROMISINGCRYSTALEXTRACTIONPARAMETERS_H
