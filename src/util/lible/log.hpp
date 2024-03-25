#pragma once

#include <string>

namespace lible
{
	namespace log
	{
		class Logger
		{
		public:

			static void writeToLog(const std::string &message);

			static bool getWriteLog()
			{
				return write_log;
			}

			static std::string getLogFName()
			{
				return log_fname;
			}

			static void setLogFName(const std::string &log_fname_in)
			{
				log_fname = log_fname_in;
			}

			static void setWriteLog(bool write_log_in)
			{
				write_log = write_log_in;
			}

		private:
			static inline bool write_log = true;
			static inline std::string log_fname = "lible.log";
		};
	}
}
