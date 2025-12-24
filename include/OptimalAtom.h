#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALATOM_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALATOM_H

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "StreamWriter.h"

#include "CrystallographicAtom.h"

#include "CoordinationConstraints.h"
#include "IonicAtomicNumber.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			class ConstrainingAtom;


			class OptimalAtom :public ChemToolkit::Crystallography::CrystallographicAtom
			{
				using size_type = unsigned short;
				using CoordinationConstraints = MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				OptimalAtom() noexcept;
				explicit OptimalAtom(const IonicAtomicNumber&) noexcept;
				OptimalAtom(const IonicAtomicNumber&, const NumericalVector&) noexcept;

				OptimalAtom(const ConstrainingAtom&) noexcept;

				virtual ~OptimalAtom() = default;

				OptimalAtom(const OptimalAtom&) = default;
				OptimalAtom(OptimalAtom&&) noexcept = default;
				OptimalAtom& operator=(const OptimalAtom&) = default;
				OptimalAtom& operator=(OptimalAtom&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const std::vector<IonicAtomicNumber>& coordinations() const noexcept;
				const std::vector<std::pair<IonicAtomicNumber, IonicAtomicNumber>>& vertexSharings() const noexcept;
				const std::vector<std::pair<IonicAtomicNumber, std::array<IonicAtomicNumber, 2>>>& edgeSharings() const noexcept;
				const std::vector<std::pair<IonicAtomicNumber, std::vector<IonicAtomicNumber>>>& faceSharings() const noexcept;

				void addCoordination(const IonicAtomicNumber&) noexcept;
				void addVertexSharings(const std::pair<IonicAtomicNumber, IonicAtomicNumber>&) noexcept;
				void addEdgeSharings(const std::pair<IonicAtomicNumber, std::array<IonicAtomicNumber, 2>>&) noexcept;
				void addFaceSharings(const std::pair<IonicAtomicNumber, std::vector<IonicAtomicNumber>>&) noexcept;

				void clearCoordinations() noexcept;
				void clearVertexSharings() noexcept;
				void clearEdgeSharings() noexcept;
				void clearFaceSharings() noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				size_type getCoordinationNumber() const noexcept;
				size_type countVertexSharings() const noexcept;
				size_type countEdgeSharings() const noexcept;
				size_type countFaceSharings() const noexcept;

				std::string getAbridgedFingerprint() const;
				std::string getFingerprint() const;
				std::string getDetailedFingerprint() const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void writeCoordinations(System::IO::StreamWriter&) const;
				void writeVertexSharings(System::IO::StreamWriter&) const;
				void writeEdgeSharings(System::IO::StreamWriter&) const;
				void writeFaceSharings(System::IO::StreamWriter&) const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				std::vector<IonicAtomicNumber> _coordinations;
				std::vector<std::pair<IonicAtomicNumber, IonicAtomicNumber>> _vertexSharings;
				std::vector<std::pair<IonicAtomicNumber, std::array<IonicAtomicNumber, 2>>> _edgeSharings;
				std::vector<std::pair<IonicAtomicNumber, std::vector<IonicAtomicNumber>>> _faceSharings;
			};



			inline bool operator<(const OptimalAtom& oaa, const OptimalAtom& oab)
			{
				const ChemToolkit::Crystallography::CrystallographicAtom& caa = dynamic_cast<const ChemToolkit::Crystallography::CrystallographicAtom&>(oaa);
				const ChemToolkit::Crystallography::CrystallographicAtom& cab = dynamic_cast<const ChemToolkit::Crystallography::CrystallographicAtom&>(oab);


				if (caa < cab)
					return true;

				else if (cab < caa)
					return false;

				else
				{
					if (oaa.getCoordinationNumber() < oab.getCoordinationNumber())
						return true;
					else if (oab.getCoordinationNumber() < oaa.getCoordinationNumber())
						return false;

					else
					{
						if (oaa.countVertexSharings() < oab.countVertexSharings())
							return true;
						else if (oab.countVertexSharings() < oaa.countVertexSharings())
							return false;

						else
						{
							if (oaa.countEdgeSharings() < oab.countEdgeSharings())
								return true;
							else if (oab.countEdgeSharings() < oaa.countEdgeSharings())
								return false;

							else
							{
								if (oaa.countFaceSharings() < oab.countFaceSharings())
									return true;
								else
									return false;
							}
						}
					}
				}
			}

			inline bool operator>(const OptimalAtom& oaa, const OptimalAtom& oab)
			{
				return (oab < oaa);
			}

			inline bool operator==(const OptimalAtom& oaa, const OptimalAtom& oab)
			{
				const ChemToolkit::Crystallography::CrystallographicAtom& caa = dynamic_cast<const ChemToolkit::Crystallography::CrystallographicAtom&>(oaa);
				const ChemToolkit::Crystallography::CrystallographicAtom& cab = dynamic_cast<const ChemToolkit::Crystallography::CrystallographicAtom&>(oab);


					if (caa == cab)
					{
						if (oaa.coordinations() == oab.coordinations())
						{
							if (oaa.vertexSharings() == oab.vertexSharings())
							{
								if (oaa.edgeSharings() == oab.edgeSharings())
								{
									if (oaa.faceSharings() == oab.faceSharings())
										return true;
									else
										throw System::ExceptionServices::InvalidOperationException{ "MathematicalCrystalChemistry::CrystalModel::Components::operator==", "Could not compare optimal atoms since number of face sharings is different." };
								}

								else
									throw System::ExceptionServices::InvalidOperationException{ "MathematicalCrystalChemistry::CrystalModel::Components::operator==", "Could not compare optimal atoms since number of edge sharings is different." };
							}

							else
								throw System::ExceptionServices::InvalidOperationException{ "MathematicalCrystalChemistry::CrystalModel::Components::operator==", "Could not compare optimal atoms since number of vertex sharings is different." };
						}

						else
							throw System::ExceptionServices::InvalidOperationException{ "MathematicalCrystalChemistry::CrystalModel::Components::operator==", "Could not compare optimal atoms since coordination is different." };
					}

					else
						return false;
			}

			inline bool operator!=(const OptimalAtom& oaa, const OptimalAtom& oab)
			{
				return (oaa != oab);
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

inline const std::vector<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber>& MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::coordinations() const noexcept
{
	return _coordinations;
}

inline const std::vector<std::pair<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber, MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber>>& MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::vertexSharings() const noexcept
{
	return _vertexSharings;
}

inline const std::vector<std::pair<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber, std::array<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber, 2>>>& MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::edgeSharings() const noexcept
{
	return _edgeSharings;
}

inline const std::vector<std::pair<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber, std::vector<MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::IonicAtomicNumber>>>& MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::faceSharings() const noexcept
{
	return _faceSharings;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::clearCoordinations() noexcept
{
	_coordinations.clear();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::clearVertexSharings() noexcept
{
	_vertexSharings.clear();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::clearEdgeSharings() noexcept
{
	_edgeSharings.clear();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::clearFaceSharings() noexcept
{
	_faceSharings.clear();
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

inline MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::getCoordinationNumber() const noexcept
{
	return static_cast<size_type>(_coordinations.size());
}

inline MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::countVertexSharings() const noexcept
{
	return static_cast<size_type>(_vertexSharings.size());
}

inline MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::countEdgeSharings() const noexcept
{
	return static_cast<size_type>(_edgeSharings.size());
}

inline MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom::countFaceSharings() const noexcept
{
	return static_cast<size_type>(_faceSharings.size());
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALATOM_H
