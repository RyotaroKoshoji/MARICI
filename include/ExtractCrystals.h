#ifndef MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_INTERNAL_ANALYSISCRYSTALS_H
#define MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_INTERNAL_ANALYSISCRYSTALS_H

#include <string>
#include <filesystem>
#include <queue>
#include <mutex>

#include "CrystalOptimalityAnalyzer.h"
#include "IsotypicCrystalExtractor.h"
#include "PromisingCrystalExtractor.h"

#include "CrystalExtractionTask.h"


namespace MathematicalCrystalChemistry
{
	namespace Extraction
	{
		namespace Internal
		{
			class ExtractCrystals
			{
				using size_type = std::size_t;
				using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;

				using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
				using SpaceGroupNumber = ChemToolkit::Crystallography::Symmetry::SpaceGroupNumber;
				using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;

				using OptimalAtom = MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom;
				using OptimalCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure;

				using CrystalOptimalityAnalyzer = MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer;
				using IsotypicCrystalExtractor = MathematicalCrystalChemistry::Analysis::IsotypicCrystalExtractor;
				using PromisingCrystalExtractor = MathematicalCrystalChemistry::Analysis::PromisingCrystalExtractor;


// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				explicit ExtractCrystals(const size_type threadIndex) noexcept;
				virtual ~ExtractCrystals() = default;

				ExtractCrystals(const ExtractCrystals&) = default;
				ExtractCrystals(ExtractCrystals&&) noexcept = default;
				ExtractCrystals& operator=(const ExtractCrystals&) = default;
				ExtractCrystals& operator=(ExtractCrystals&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				void setCrystalExtractionTask(const CrystalExtractionTask&) noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void operator()() const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				std::queue<std::filesystem::path> getCifFilePaths() const;
				void eraseFilePathsInitially(std::queue<std::filesystem::path>&) const;
				void eraseFilePathsForParallel(std::queue<std::filesystem::path>&) const;

				void outputOptimalCrystalStructure(std::filesystem::path inputCifFilePath, const OptimalCrystalStructure&) const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private utility

				std::filesystem::path getSpaceGroupDirectoryPath(const SpaceGroupNumber, const ChemicalComposition&) const;
				std::pair<bool, std::filesystem::path> getOutputDirectoryPath(const std::string& outputFingerprint, const std::filesystem::path& spaceGroupPath) const;

				void outputFeasibleCrystallographicData(const OptimalCrystalStructure& conventionalOptimalStructure, const std::filesystem::path& outputDirectoryPath) const;

			// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				size_type _threadRank;
				std::filesystem::path _inputCrystalsDirectoryPath;
				std::filesystem::path _outputCrystalsDirectoryPath;

				double _spaceGroupPrecision;
				double _feasibleErrorRate;
				GeometricalConstraintParameters _geometricalConstraintParameters;

				CrystalOptimalityAnalyzer _crystalOptimalityAnalyzer;
				IsotypicCrystalExtractor _isotypicCrystalExtractor;
				PromisingCrystalExtractor _promisingCrystalExtractor;


				static std::string s_fingerprintFilename;
				static std::mutex s_outputMutex;
				static double s_defaultFeasibleErrorRate;
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

inline void MathematicalCrystalChemistry::Extraction::Internal::ExtractCrystals::setCrystalExtractionTask(const CrystalExtractionTask& task) noexcept
{
	_inputCrystalsDirectoryPath = task.inputCrystalsDirectoryPath();
	_outputCrystalsDirectoryPath = task.outputCrystalsDirectoryPath();

	_spaceGroupPrecision = task.spaceGroupPrecision();
	_feasibleErrorRate = task.crystalOptimalityAnalysisParameters().preciseStructuralOptimizationParameters().feasibleGeometricalConstraintErrorRate();
	_geometricalConstraintParameters = task.geometricalConstraintParameters();

	_crystalOptimalityAnalyzer.setOptimalityAnalysisParameters(task.crystalOptimalityAnalysisParameters());
	_isotypicCrystalExtractor.setIsotypicCrystalExtractionParameters(task.isotypicCrystalExtractionParameters());
	_promisingCrystalExtractor.setPromisingCrystalExtractionParameters(task.promisingCrystalExtractionParameters());
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_EXTRACTION_INTERNAL_ANALYSISCRYSTALS_H
