class PipelineBase:

    def __init__(self):
        self.flags = {}  # bitmask flag

    def set_flag(self, name):
        """Set a flag using bitmask"""
        #self.flags |= value
        self.flags[name] = True 

    def clear_flag(self, value):
        """Remove a flag"""
        self.flags &= ~value

    def check_flag(self, value):
        """Check if flag is set"""
        return (self.flags & value) != 0


    def log_issue(self, func_name, level, message):
        with open("logging_pymorph.log", "a") as f:
            name = getattr(self, "NAME", "UNKNOWN")
            f.write(f"{name} | {func_name} | {level} | {message}\n")
