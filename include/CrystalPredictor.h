#ifndef MATHEMATICALCRYSTALCHEMISTRY_PREDICTION_CRYSTALPREDICTOR_H
#define MATHEMATICALCRYSTALCHEMISTRY_PREDICTION_CRYSTALPREDICTOR_H

#include <mutex>

#include "ConstrainingAtomicSpecies.h"

#include "RandomStructureGenerator.h"
#include "CrystalProductionReporter.h"
#include "CrystalDesigner.h"

#include "CrystalPredictionTask.h"


namespace MathematicalCrystalChemistry
{
	namespace Prediction
	{
		class CrystalPredictor
		{
			using size_type = std::size_t;
			using ConstrainingAtomicSpecies = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies;
			using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<ConstrainingAtomicSpecies>;

			using CrystalDesigner = MathematicalCrystalChemistry::Design::CrystalDesigner;
			using RandomStructureGenerator = MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerator;
			using CrystalProductionReporter = MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReporter;

			using CrystalDesignParameters = MathematicalCrystalChemistry::Design::CrystalDesignParameters;
			using CrystalProductionReportParameters = MathematicalCrystalChemistry::Design::Diagnostics::CrystalProductionReportParameters;


			class ProduceCrystals
			{
			public:
				explicit ProduceCrystals(const size_type threadRank) noexcept;
				virtual ~ProduceCrystals() = default;

				ProduceCrystals(const ProduceCrystals&) = default;
				ProduceCrystals(ProduceCrystals&&) noexcept = default;
				ProduceCrystals& operator=(const ProduceCrystals&) = default;
				ProduceCrystals& operator=(ProduceCrystals&&) noexcept = default;


				void setStdFilePath() noexcept;
				void setStdFilePath(const std::filesystem::path&);
				void setCrystalDesignParameters(const CrystalDesignParameters&);
				void setCrystalProductionReportParameters(const CrystalProductionReportParameters&, const std::filesystem::path& mpiCrystalProductionDirectoryPath);

				void operator()();

				static void initializeStructureProducing();
				static void setMaxStructureProducing(const size_type);
				static void setChemicalComposition(const ChemicalComposition&);


			private:
				bool shouldDesign() const;
				void reportCeaselessGeneration() const;

				std::string getProductionName() const;
				ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> toBasicChemicalComposition(const ChemicalComposition&) const;


			private:
				std::filesystem::path _stdFilePath;

				size_type _threadRank;
				RandomStructureGenerator _randomStructureGenerator;
				CrystalDesigner _crystalDesigner;
				CrystalProductionReporter _crystalProductionReporter;

				mutable size_type m_sampleIndex;

				static ChemicalComposition s_chemicalComposition;
				static size_type s_structureDesigning;
				static size_type s_maxStructureDesigning;
				static std::mutex s_productionMutex;
			};

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalPredictor() noexcept;
			explicit CrystalPredictor(const CrystalPredictionTask&);

			virtual ~CrystalPredictor() = default;

			CrystalPredictor(const CrystalPredictor&) = default;
			CrystalPredictor(CrystalPredictor&&) noexcept = default;
			CrystalPredictor& operator=(const CrystalPredictor&) = default;
			CrystalPredictor& operator=(CrystalPredictor&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			bool isStdoutEnabled() const noexcept;

			void enableStdout() noexcept;
			void disableStdout() noexcept;
			void setCrystalPredictionTask(const CrystalPredictionTask&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			void execute() const;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Private methods

		private:
			void updateStdFilePath(const std::filesystem::path& crystalProductionDirectoryPath) const;

			void createMpiProductionDirectoryPaths(const std::filesystem::path& crystalProductionDirectoryPath) const;
			void createLogFiles(const std::filesystem::path& crystalProductionDirectoryPath) const;
			void initializeCrystalProducers(const std::filesystem::path& crystalProductionDirectoryPath) const;

			void validateProductionPaths(const std::filesystem::path& crystalProductionDirectoryPath) const;

			void reportAllJobFinalization() const;
			void reportJobFinalization(const ChemicalComposition&, const size_type) const;
			size_type getMaxCrystalProducing(const size_type requestedCrystalProducing) const;

		// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			bool _isStdoutEnabled;
			CrystalPredictionTask _crystalPredictionTask;

			mutable std::filesystem::path m_stdFilePath;
			mutable std::vector<ProduceCrystals> m_crystalProducers;
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

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setStdFilePath() noexcept
{
	_stdFilePath.clear();
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setStdFilePath(const std::filesystem::path& filePath)
{
	if (filePath.empty())
		throw System::ExceptionServices::ArgumentOutOfRangeException(typeid("this"), "setStdFilePath", "Argument file path is empty.");
	else
		_stdFilePath = filePath;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setCrystalDesignParameters(const CrystalDesignParameters& parameters)
{
	_randomStructureGenerator.setRandomStructureGenerationParameters(parameters.initialStructureGenerationParameters().randomStructureGenerationParameters());
	_crystalDesigner.setCrystalDesignParameters(parameters);
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setCrystalProductionReportParameters(const CrystalProductionReportParameters& parameters, const std::filesystem::path& mpiCrystalProductionDirectoryPath)
{
	_crystalProductionReporter.setCrystalProductionReportParameters(parameters);
	_crystalProductionReporter.crystalDesignRecorder().setMpiCrystalProductionDirectoryPath(mpiCrystalProductionDirectoryPath);
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::initializeStructureProducing()
{
	std::lock_guard<std::mutex> guard{ s_productionMutex };
	s_structureDesigning = 0;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setMaxStructureProducing(const size_type num)
{
	std::lock_guard<std::mutex> guard{ s_productionMutex };
	s_maxStructureDesigning = num;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::setChemicalComposition(const ChemicalComposition& composition)
{
	std::lock_guard<std::mutex> guard{ s_productionMutex };
	s_chemicalComposition = composition;
}


inline bool MathematicalCrystalChemistry::Prediction::CrystalPredictor::isStdoutEnabled() const noexcept
{
	return _isStdoutEnabled;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::enableStdout() noexcept
{
	_isStdoutEnabled = true;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::disableStdout() noexcept
{
	_isStdoutEnabled = false;
}

inline void MathematicalCrystalChemistry::Prediction::CrystalPredictor::setCrystalPredictionTask(const CrystalPredictionTask& task) noexcept
{
	_crystalPredictionTask = task;
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

inline bool MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::shouldDesign() const
{
	std::lock_guard<std::mutex> guard{ s_productionMutex };

	if (s_structureDesigning < s_maxStructureDesigning)
	{
		m_sampleIndex = ++s_structureDesigning;
		return true;
	}

	else
		return false;
}

inline ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> MathematicalCrystalChemistry::Prediction::CrystalPredictor::ProduceCrystals::toBasicChemicalComposition(const ChemicalComposition& chemicalComposition) const
{
	ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> basicComposition;
	{
		for (const auto& speciesAndCount : chemicalComposition)
			basicComposition.add(speciesAndCount.first.ionicAtomicNumber().atomicNumber(), speciesAndCount.second);
	}

	return basicComposition;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_PREDICTION_CRYSTALPREDICTIONTOR_H
