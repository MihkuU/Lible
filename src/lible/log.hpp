#pragma once

#include <string>
#include <fstream>

namespace lible
{
	namespace log
	{
		class Logger
		{
		public:

			Logger();

			~Logger();

			void operator<<(const std::string &message);

		private:
			std::string log_fname;

			std::ofstream log_file;
		};

		inline Logger logger{};
	}
}
