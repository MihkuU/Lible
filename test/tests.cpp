#include <tests.hpp>

#include <lible/ints/ints.hpp>

namespace ltests = lible::tests;
namespace lints = lible::ints;

bool ltests::numCartesians()
{
	int correct_answer = 21;

	int l = 5;
	int num_cartesians = lints::numCartesians(l);

	if (num_cartesians == correct_answer)
		return true;
	else 
		return false;
}

bool ltests::numSphericals()
{
	int correct_answer = 13;

	int l = 6;
	int num_sphericals = lints::numSphericals(l);

	if (num_sphericals == correct_answer)
		return true;
	else
		return false;
}
