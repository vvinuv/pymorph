from functools import wraps
from pathlib import Path
import warnings
import math
import numpy as np


class PipelineBase:

    def __init__(self):
        self.flags = {}  # bitmask flag

    def set_flag(self, step, level, msg):
        """Set a flag using bitmask"""
        self.flags[step] = {"level": level, "msg": msg}
        #self.flags |= value
        #self.flags[name] = True

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
            "error": "CRITICAL",
            "issue": "LOG_NOT_LIST"
        }


    if not isinstance(lst, (list, tuple)):
        print('NOT_A_LIST')
        return False, {
            "error": "CRITICAL",
            "issue": "NOT_A_LIST"
        }

    if any(isinstance(x, (list, tuple)) for x in lst):
        print('NESTED_LIST')
        return False, {
            "error": "CRITICAL",
            "issue": "NESTED_LIST"
        }

    if len(lst) != expected_len:
        print('INVALID_LENGTH')
        return False, {
            "error": "CRITICAL",
            "issue": "INVALID_LENGTH"
        }

    #print('What the fuck 3 True')
    return True, {}


def check_params_and_errors(params, errors, expected_length):

    #print('params', params)
    #print('expected_length', expected_length)
    #  ---- CHECK PARAMS (CRITICAL) ----
    #if (not params or
    #    any(p is None or (isinstance(p, float) and np.isnan(p)) or np.isinf(p) for p in params)):
    if len(params) != expected_length:
        #print('What the fuck INVALID_PARAMETERS') 
        return False, (), {
            "error": "CRITICAL",
            "issue": "INVALID_PARAMETERS"
        }

    #  ---- FIX ERRORS ----
    # enforce same length always
    if (not errors) or (len(errors) != len(params)):
        errors = [9999] * len(params)

        return True, (params, errors), {
                "error": "Warning: ERRORS_REPLACED",
                "issue": "ERRORS_EMPTY_OR_SIZE_MISMATCH"
                }

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
                "error": "Warning: ERROR_NAN_VALUES to 9999",
                "issue": "NAN_VALUES_REPLACED_9999"
                }

    return True, (params, errors), {"error": "", "issue": ""}


def check_image_size(value):
    # EXAMPLE RULE: THE VALUE SHOULD GREATER THAN 50
    if value < 30:
        return False, {
            "error": "CRITICAL",
            "issue": "SMALL_IMAGE_SIZE"
            }

    return True, {}



def check_file(fname):

    if not Path(fname).exists():
        return False, fname, {
            "error": "CRITICAL",
            "issue": "GALFIT_FAILED"
        }

    return True, fname, {}



def catch_pipeline_issues(critical=None, file_checker=False, 
                          expected_len=None, value_checker=None):

    def decorator(func):

        @wraps(func)
        def wrapper(self, *args, **kwargs):

            # reset
            #self.flags = {}

            # RUN FUNCTION
            result = func(self, *args, **kwargs)
            #print(f'result  {result} {func.__name__}')
            #if len(result) != expected_len:

            # GENERIC OUTPUT CHECK. EARLIER IT IS THE FUNCTION INTO THIS
            # DECORDOR. NOW IT IS USING THE ACTUAL FUNCTION WHICH GOT 
            # THREE OUTPUTS

            # SIZE CHECK
            if isinstance(result, int) and result < 30:
                #HERE RESULT IS THE SIZE OF THE IMAGE

                status, info = check_image_size(result)

                if not status:
                    info = {
                        "function": func.__name__,
                        **info
                    }
                    self.set_flag("SMALL_IMAGE_SIZE", 0)
                    raise PipelineCriticalError({
                                            "error": info["error"],
                                            "issue": "SMALL_IMAGE_SIZE",
                                        })
                return result

            #CHECK WHETHER FIT.LOG EXISTS
            if file_checker == True: 
                #HERE RESULT IS THE FILE NAME
                status, fname, info = check_file(result)
                if not status:
                    info = {
                            "function": func.__name__,
                            **info
                            }
                    f = fname.split('.')[0].upper()
                    self.set_flag(f"GALFIT_FAILED", 1)
                    raise PipelineCriticalError({
                                            "error": info["error"],
                                            "issue": "GALFIT_FAILED"
                                        })
                return info 

            
            #CHECK WHETHER THE EXPECTED LENGTH IS SATIFIED. 
            #IT IS NOT IMPORTANT NOW. USE CHECK_PARAMS_AND_ERRORS FUNCTION 
            #INSTEAD OF THIS. HOWEVER, IF WE HAVE A LIST WHICH HAS DIFFERENT
            #LENGTH INSTEAD OF AN EXPECTED LENGTH OR NAN OR INF OR NONE WILL 
            #BE TAKEN CARE OF. ANYWAY IT IS NOT IMPORTANT NOW.
            if expected_len is not None and 0:
                status, info = find_len(result, expected_len)

                if not status:
                    info = {
                        "function": func.__name__,
                        **info
                    }
                    raise PipelineCriticalError({
                        "error": info["error"],
                        "issue": "PARAMS_LEN_LESS_DEFAULT"
                        })


            # PARAMS + ERRORS CHECK (SEPARATE FUNCTION INSTEAD OF THE ABOVE)
            if isinstance(result, tuple) and len(result) == 2:
                #print("PARAMS + ERRORS")
                params, errors = result
                #print(result)
                status, result, info = check_params_and_errors(params,
                                                               errors,
                                                               expected_len)

                if not status:
                    info = {
                        "function": func.__name__,
                        **info
                    }
                    self.pipeline.set_flag("PARAMS_NO_FLOAT", 2)

                    raise PipelineCriticalError({
                        "error": info["error"],
                        "issue": "PARAMS_NO_FLOAT"
                        })
                
                if status:
                    info = {
                        "function": func.__name__,
                        **info
                    }
                    if info["issue"] != "":  
                        self.pipeline.set_flag(info["issue"], 5)


                return result[1]

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
#    def check_image_size(self, value):
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

class MyClass(PipelineBase):

    def __init__(self):
        self.flags = {}

    def set_flag(self, name, flag):
        if name not in self.flags:
            self.flags[name] = flag   # create list


    def compute_list(self, lst):
        return lst   # you can change this for testing


    # uses ONLY value check
    @catch_pipeline_issues()
    def check_image_size(self, value):
        return value

    @catch_pipeline_issues(expected_len=7)
    def check_params_and_errors(self, params, errors):
        return params, errors

    @catch_pipeline_issues(file_checker=True)
    def check_file(self, filename):
        return filename

    def parse(self, lst):
        
        #try:

        self.set_flag('BOX', 3)
        self.set_flag('BAR', 4)
        print('flags', self.flags)

        Ifile = self.check_file('t.txt')

        size = self.check_image_size(30)
        print('flags', self.flags)

        #print('size', size)
        params = self.compute_list(lst)
       
        #print('Inside1', params)
        #params = [1, 2, 3, 4, 5, 6, np.inf]
        
        param_err = self.check_params_and_errors(params, [1, 2]) 
        print('Inside', param_err)
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
    fla = PipelineBase() 
    for i in range(2):
        #lst.append(i)
        obj = MyClass()
        print(i)
        try:
            obj.parse(lst)
            print('obj.flags', obj.flags)
        except PipelineCriticalError as e:
            print(obj.flags)
            print("Caught:", e)
            continue   # ✅ go to next iteration:


#obj.compute_list(3)
#obj.check_image_size(-2)
#print(obj.critical)
#print(obj.flags)




#lst = [2, 3]
#obj = MyClass_prev()
#obj.parse(lst)
#print(obj.critical)
#print(obj.flags)



