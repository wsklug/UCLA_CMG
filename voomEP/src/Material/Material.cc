//-*-C++-*-
#include "Material.h"
#include "CompNeoHookean.h"
#include "CompMyoCardium.h"
#include "NeoHookean.h"

namespace voom{
  //! Gateway to material class
  Material* Material::New(Real data[], materialType type) {
    Material* myMaterial;
    switch(type) {
    case 0:
      myMaterial = new CompNeoHookean(data);
      break;
    case 1:
      myMaterial = new CompMyoCardium(data);
      break;
    case 2:
      myMaterial = new NeoHookean(data);
      break;
    default:
      {
	std::cout << "Unknown Material Type specified\n";
	exit(1);
      }
    }
    return myMaterial;
  }

} // namespace

