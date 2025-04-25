#include <format>
#include <stdexcept>

#include <tests.hpp>

int main(int argc, char** argv)
{
    if (argc != 2)
	{
		std::string msg = std::format("Invalid number of arguments: {}. Only one argument must be given!", (argc - 1));
      	throw std::runtime_error(msg);
	}

    std::string test_name = argv[1];	

	bool success = false;
	if (test_name == "numCartesians")
		success = lible::tests::numCartesians();
	else if (test_name == "numSphericals")
		success = lible::tests::numSphericals();
	else 
		throw std::runtime_error(std::format("Invalid test name specified: {}.", test_name));

    if (success)
        return 0;
    else
        return 1;
}
