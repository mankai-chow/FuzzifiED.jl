"""
    SilentStd :: bool = false 

a flag to determine whether logs of the FuzzifiED functions should be turned off. False by default. If you want to evaluate without log, put `FuzzifiED.SilentStd = true`
"""
SilentStd = false
macro ctrlstd(ex)
    quote
        @show SilentStd
        if SilentStd
            # Save the current stdout and stderr
            original_stdout = stdout
            original_stderr = stderr
            
            # Redirect stdout and stderr to /dev/null
            open("/dev/null", "w") do devnull
                redirect_stdout(devnull)
                redirect_stderr(devnull)
                
                try
                    # Evaluate the expression silently
                    $(esc(ex))
                finally
                    # Restore the original stdout and stderr
                    redirect_stdout(original_stdout)
                    redirect_stderr(original_stderr)
                end
            end
        else
            # Evaluate the expression normally
            $(esc(ex))
        end
    end
end