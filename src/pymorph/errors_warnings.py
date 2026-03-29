from functools import wraps
from pathlib import Path
import warnings
import math
import numpy as np


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




class PipelineCriticalError(Exception):
    #pass
    def __init__(self, info):
        self.info = info
        super().__init__(str(info))


def find_len(lst, expected_len):
    #print('What the fuck 2', lst, expected_len)

    if not lst:
        print('EMPTY')
        return False, {
            "issue": "EMPTY"
        }, "CRITICAL"

    if not isinstance(lst, (list, tuple)):
        print('NOT_A_LIST')
        return False, {
            "issue": "NOT_A_LIST",
            "actual_type": type(lst).__name__
        }, "CRITICAL"

    if any(isinstance(x, (list, tuple)) for x in lst):
        print('NESTED_LIST')
        return False, {
            "issue": "NESTED_LIST",
            "actual_type": type(lst).__name__
        }, "CRITICAL"


    if len(lst) != expected_len:
        print('INVALID_LENGTH')
        return False, {
            "issue": "INVALID_LENGTH",
            "expected": expected_len,
            "actual": len(lst)
        }, "CRITICAL"

    #print('What the fuck 3 True')
    return True, {}, None


def check_params_and_errors(params, errors):

    #print('params', params)
    #  ---- CHECK PARAMS (CRITICAL) ----
    if (not params or
        any(p is None or (isinstance(p, float) and np.isnan(p)) or np.isinf(p) for p in params)):

        #print('What the fuck INVALID_PARAMETERS') 
        return False, (params, errors), {
            "stage": "params",
            "issue": "INVALID_PARAMETERS",
            "params": params
        }, "CRITICAL"

    #  ---- FIX ERRORS ----
    # enforce same length always
    if (not errors) or (len(errors) != len(params)):
        errors = [9999] * len(params)

        return True, (params, errors), {
            "stage": "errors",
            "issue": "ERRORS_REPLACED",
            "reason": "EMPTY_OR_SIZE_MISMATCH",
            "expected_len": len(params)
        }, "WARNING"

    # replace NaNs in errors
    new_errors = []
    replaced = False

    for e in errors:
        if e is None or (isinstance(e, float) and math.isnan(e)):
            new_errors.append(9999)
            replaced = True
        else:
            new_errors.append(e)

    if replaced:
        errors = new_errors

        return True, (params, errors), {
            "stage": "errors",
            "issue": "ERRORS_REPLACED",
            "reason": "NAN_VALUES"
        }, "WARNING"

    return True, (params, errors), {}, None


def check_value(value):
    # EXAMPLE RULE: THE VALUE SHOULD GREATER THAN 50
    if value < 50:
        return False, {
            "issue": "IMAGE_SIZE_SMALL",
            "value": value
        }, "CRITICAL"

    return True, {}, None



def check_file(filename, flag):

    if not Path(filename).exists():
        return False, (filename, flag), {
            "issue": "FILE_NOT_FOUND",
            "value": filename,
            "flag": flag,
        }, "CRITICAL"

    return True, (filename, flag), {}, None



def catch_pipeline_issues(critical=None, file_checker=False, 
                          expected_len=None, value_checker=None):

    def decorator(func):

        @wraps(func)
        def wrapper(self, *args, **kwargs):

            # reset
            self.flags = {}
            self.critical = False

            # RUN FUNCTION
            result = func(self, *args, **kwargs)

            # GENERIC OUTPUT CHECK. EARLIER IT IS THE FUNCTION INTO THIS
            # DECORDOR. NOW IT IS USING THE ACTUAL FUNCTION WHICH GOT 
            # THREE OUTPUTS

            if file_checker == True and isinstance(result, tuple) and len(result) == 2: 
                filename = result[0]
                flag = result[1]
                status, (filename, flag), info, severity = check_file(filename, flag)
                if not status:
                    self.flags = {
                        "function": func.__name__,
                        "stage": "input",
                        **info
                    }
                    self.critical = self.flags["flag"]
                    raise PipelineCriticalError({
                                            "reason": self.flags["value"],
                                            "issue": "NOT_FOUND",
                                            "flag": self.flags["flag"]
                                        })
                    return None
                return info 

            # SIZE CHECK
            if isinstance(result, int) and result < 50:

                status, info, severity = check_value(result)

                if not status:
                    self.flags = {
                        "function": func.__name__,
                        "stage": "input",
                        **info
                    }
                    self.critical = True
                    raise PipelineCriticalError("FIRST_ERROR")
                    return None
                return info["value"]

            
            if expected_len is not None:
                status, info, severity = find_len(result, expected_len)

                if not status:
                    self.flags = {
                        "function": func.__name__,
                        "stage": "output",
                        **info
                    }
                    self.critical = True
                    raise PipelineCriticalError("FIRST_ERROR")
                    return None


            # PARAMS + ERRORS CHECK (separate function)
            if isinstance(result, tuple) and len(result) == 2:
                print("PARAMS + ERRORS")
                params, errors = result
                #print(result)
                status, (params, errors), info, severity = check_params_and_errors(params, errors)

                if info:
                    self.flags = {
                        "function": func.__name__,
                        **info
                    }

                if not status:
                    self.critical = True
                    raise PipelineCriticalError("FIRST_ERROR")
                    return None

                result = (params, errors)
                
            
            return result

        return wrapper
    return decorator




#class MyClass_prev:
#
#    def __init__(self):
#        self.flags = {}
#        self.critical = False
#
#    # uses BOTH checks
#    @catch_pipeline_issues(expected_len=5)
#    def compute_list(self, lst):
#        return lst
#
#    # uses ONLY value check
#    @catch_pipeline_issues(value_checker=2)
#    def compute_value(self, value):
#        return value * 2
#
#    def parse(self, value):
#        try:
#            result = self.compute_list(value)
#            return result
#
#        except PipelineCriticalError as e:
#            print("Handled in parse:", self.flags)
#            return None
#            #print("[CRITICAL ERROR CAUGHT IN PARSE]")
#            #print("Details:", e)
#            #return None
#

class MyClass:

    def __init__(self):
        self.flags = {}
        self.critical = False

    @catch_pipeline_issues(expected_len=7)#, value_checker=check_value)
    def compute_list(self, lst):
        return lst   # you can change this for testing


    # uses ONLY value check
    @catch_pipeline_issues()
    def compute_value(self, value):
        return value

    @catch_pipeline_issues()
    def check_params_and_errors(self, params, errors):
        return params, errors

    @catch_pipeline_issues(file_checker=True)
    def check_file(self, filename, flag):
        return filename

    def parse(self, lst):
        
        #try:

        Ifile = self.check_file('t.txt', 1)

        size = self.compute_value(50)

        #print('size', size)
        params = self.compute_list(lst)
        
        #print('Inside1', params)
        #params = [1, 2, 3, 4, 5, 6, np.inf]
        
        param_err = self.check_params_and_errors(params, [1, 2]) 
        #print('Inside', param_err)
        #return result


        #except PipelineCriticalError as e:
        #    print("[CRITICAL ERROR CAUGHT IN PARSE]")
        #    print("Details:", e)
        #    return None

if __name__=='__main__':

    lst = [1, 2, 3, 4, 5, 6, 7]

##obj = MyClass()
##result =obj.parse(lst)
##print(result)
#lst = []
    for i in range(1):
        #lst.append(i)
        obj = MyClass()
        print(i)
        try:
            obj.parse(lst)
        except PipelineCriticalError as e:
            print("Caught:", e)
            continue   # ✅ go to next iteration:


#obj.compute_list(3)
#obj.compute_value(-2)
#print(obj.critical)
#print(obj.flags)




#lst = [2, 3]
#obj = MyClass_prev()
#obj.parse(lst)
#print(obj.critical)
#print(obj.flags)



