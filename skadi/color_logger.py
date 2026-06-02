import sys
import logging
import traceback


# ANSI escape sequences for colored output
COLORS = {
    logging.DEBUG: "\x1b[38;20m",  # Grey
    logging.INFO: "\x1b[37;20m",  # Green
    logging.WARNING: "\x1b[33;20m",  # Yellow
    logging.ERROR: "\x1b[31;20m",  # Red
    logging.CRITICAL: "\x1b[31;1m",  # Bold Red
}
RESET = "\x1b[0m"

PREFIX = "[skadi] "


class ColorFormatter(logging.Formatter):
    """Custom formatter to add colors to log levels."""

    def __init__(
        self,
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)",
    ):
        super().__init__(fmt)  # Ensure parent class initializes properly
        self.FORMATS = {
            level: PREFIX + color + fmt + RESET for level, color in COLORS.items()
        }

    def format(self, record):
        """Format log messages with the corresponding color."""
        log_fmt = self.FORMATS.get(
            record.levelno, self._fmt
        )  # Use default if not found
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


# Define logging format
LOGGING_FORMAT = "%(asctime)s | %(levelname)s: %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M"

# Create console logging handler
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.INFO)
consoleHandler.setFormatter(ColorFormatter(LOGGING_FORMAT))

# Configure the root logger
logging.basicConfig(
    level=logging.DEBUG,
    datefmt=DATE_FORMAT,
    format=LOGGING_FORMAT,
    handlers=[consoleHandler],
)

logging.captureWarnings(True)  # Capture warnings as log messages
# Root logger
logger = logging.getLogger()


def handle_exception(exc_type, exc_value, exc_traceback):
    """Global exception handler to log uncaught exceptions."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "Uncaught exception:\n%s",
        "".join(traceback.format_exception(exc_type, exc_value, exc_traceback)),
    )


# Install exception handler
sys.excepthook = handle_exception
