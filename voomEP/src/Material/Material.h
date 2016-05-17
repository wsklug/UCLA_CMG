//-*-C++-*-
// Material Class

/*!
  \brief Base class for a material objects.
*/

#ifndef _Material_h_
#define _Material_h_

#include "voom.h"
#include "VoomMath.h"

namespace voom {
  //! Specify material type to be used
  enum materialType {COMPNEOHOOKEAN, COMPMYOCARDIUM, NEOHOOKEAN};
  class Material {
  protected:
    //! Strain Energy
    Real _W;

    //! Deformation Gradient
    Tensor3D _F;

    //! 1st Piola Kirchoff Stress
    Tensor3D _P;

    //! Cauchy Stress
    Tensor3D _S;

    //! Lagrangian Modulus of the material
    Array4D _C;

    //! I1, I2, I3
    Real   _I1, _I2, _I3;

  public:
    //! Default Constructor
    Material() {;}

    //! Destructor
    ~Material() {;}

    //! Set deformation gradient
    void setDeformationGradient(const Tensor3D& F) {
      if (determinant(F) > 0. ) {
	_F = F;
      } else {
	std::cout << F << "\n";
	std::cout << determinant(F) << "\n\n";
	std::cout << "Jacobian (DefGrad) is negative Exiting" << std::endl;
	exit(1);
      }
    }

    //! Get an invariant from the model
    virtual Real getInvariant(std::string request) const = 0;

    //! Get Piola Stress
    virtual const Tensor3D& piolaStress() const { return _P; }

    //! Get Cauchy Stress
    virtual const Tensor3D& cauchyStress() const { return _S; }

    //! Get strain energy
    virtual Real strainEnergy() const { return _W; }

    //! Get Stiffness Matrix
    const Array4D& getLagrangianModulus() const { return _C; }

    //! Update state. Compute W, P and C
    virtual void updateState(bool f0, bool f1, bool f2, 
			     const Vector3D f=Vector3D(0.,0.,1.),
			     const Vector3D s=Vector3D(1.,0.,0.) ) = 0;

    /*! 
      Gateway to material class. The parameters
      \param data a Real list containing pertinent material parameters
      \param type a materialType variable which defines the type of
      material to be used
      \return A base class pointer pointing to the derived material object
    */
    static Material* New(Real data[], materialType type);
  };
} // namespace

#endif
