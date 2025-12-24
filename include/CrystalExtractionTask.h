#ifndef MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_CRYSTALEXTRACTIONTASK_H
#define MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_CRYSTALEXTRACTIONTASK_H

#include <string>
#include <filesystem>

#include "GeometricalConstraintParameters.h"

#include "CrystalOptimalityAnalysisParameters.h"
#include "IsotypicCrystalExtractionParameters.h"
#include "PromisingCrystalExtractionParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Extraction
	{
		class CrystalExtractionTask
		{
			using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;

			using CrystalOptimalityAnalysisParameters = MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters;
			using IsotypicCrystalExtractionParameters = MathematicalCrystalChemistry::Analysis::IsotypicCrystalExtractionParameters;
			using PromisingCrystalExtractionParameters = MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalExtractionTask() noexcept;
			virtual ~CrystalExtractionTask() = default;

			CrystalExtractionTask(const CrystalExtractionTask&) = default;
			CrystalExtractionTask(CrystalExtractionTask&&) noexcept = default;
			CrystalExtractionTask& operator=(const CrystalExtractionTask&) = default;
			CrystalExtractionTask& operator=(CrystalExtractionTask&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			const std::filesystem::path& inputCrystalsDirectoryPath() const noexcept;
			const std::filesystem::path& outputCrystalsDirectoryPath() const noexcept;

			double spaceGroupPrecision() const noexcept;
			const GeometricalConstraintParameters& geometricalConstraintParameters() const noexcept;

			const CrystalOptimalityAnalysisParameters& crystalOptimalityAnalysisParameters() const noexcept;
			const IsotypicCrystalExtractionParameters& isotypicCrystalExtractionParameters() const noexcept;
			const PromisingCrystalExtractionParameters& promisingCrystalExtractionParameters() const noexcept;

			static double defaultSpaceGroupPrecision() noexcept;


			void setInputDirectoryPath(const std::filesystem::path&);
			void setOutputDirectoryPath(const std::filesystem::path&);

			void setGeometricalConstraintParameters(const GeometricalConstraintParameters&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			void initialize() noexcept;
			void initialize(const System::IO::StreamReader&);

			bool hasTask() const noexcept;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Private methods

		private:
			void validateInitializedValues() const;

		// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			std::filesystem::path _inputCrystalsDirectoryPath;
			std::filesystem::path _outputCrystalsDirectoryPath;

			double _spaceGroupPrecision;
			GeometricalConstraintParameters _geometricalConstraintParameters;

			CrystalOptimalityAnalysisParameters _crystalOptimalityAnalysisParameters;
			IsotypicCrystalExtractionParameters _isotypicCrystalExtractionParameters;
			PromisingCrystalExtractionParameters _promisingCrystalExtractionParameters;

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

inline const std::filesystem::path& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::inputCrystalsDirectoryPath() const noexcept
{
	return _inputCrystalsDirectoryPath;
}

inline const std::filesystem::path& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::outputCrystalsDirectoryPath() const noexcept
{
	return _outputCrystalsDirectoryPath;
}

inline double MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::spaceGroupPrecision() const noexcept
{
	return _spaceGroupPrecision;
}

inline const MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::geometricalConstraintParameters() const noexcept
{
	return _geometricalConstraintParameters;
}

inline double MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::defaultSpaceGroupPrecision() noexcept
{
	return s_defaultSpaceGroupPrecision;
}

inline const MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::crystalOptimalityAnalysisParameters() const noexcept
{
	return _crystalOptimalityAnalysisParameters;
}

inline const MathematicalCrystalChemistry::Analysis::IsotypicCrystalExtractionParameters& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::isotypicCrystalExtractionParameters() const noexcept
{
	return _isotypicCrystalExtractionParameters;
}

inline const MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractionParameters& MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::promisingCrystalExtractionParameters() const noexcept
{
	return _promisingCrystalExtractionParameters;
}

inline void MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::setGeometricalConstraintParameters(const GeometricalConstraintParameters& parameters) noexcept
{
	_geometricalConstraintParameters = parameters;
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

inline bool MathematicalCrystalChemistry::Extraction::CrystalExtractionTask::hasTask() const noexcept
{
	if (_inputCrystalsDirectoryPath.empty())
		return false;

	if (_outputCrystalsDirectoryPath.empty())
		return false;

	if (!(_crystalOptimalityAnalysisParameters.needAnalysis()) && !(_isotypicCrystalExtractionParameters.needExtraction()) && !(_promisingCrystalExtractionParameters.needExtraction()))
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


#endif // !MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_CRYSTALEXTRACTIONTASK_H
