#! /usr/bin/env ruby
# Wrapper for surrogate algorithm call.

# Example output...
#
# number, dummy_non_parameter_1, dummy_non_parameter_2, dummy_non_parameter_3, dummy_non_parameter_4, alpha, ps, rho, wp
# 0, 50.0, 1, 1, 5.0, 1.3, 0.05, 0.8, 0.01

require 'tempfile'

def main(argv)
	outfilename_prefix = 'tmp_';	
	params = {}
	
	i = 5
	while i <= argv.length-2
		param = argv[i].sub(/^-/,"")
		params[param] = argv[i+1]
		i += 2
	end
	
	timeout = 300
	res = nil
    
    output_file = Tempfile.open(outfilename_prefix, './')
    begin
		output_file.print "number, dummy_non_parameter_1, dummy_non_parameter_2, dummy_non_parameter_3, dummy_non_parameter_4, "
		output_file.puts params.keys.sort.join(", ")
		values = []
		for key in params.keys.sort
			values << params[key]
		end
		output_file.print "0, 0, 0, 0, 0, "
		output_file.puts values.join(", ")
	
	    cmd = "./EPM_predict_saps #{output_file.path} #{argv[0]}"
	    puts cmd
	    res = nil	    
    ensure
        output_file.close
    end
    
    # Call EPM_predict_saps
    File.popen(cmd) {|file|
	    while line=file.gets
		    res = line.chomp;
	    end
    }    
    output_file.unlink
    	
	solved = "SAT";
	measured_runtime = res;
	obj = 0;
	gap = 0;
	seed = argv[4]
	puts "Result for ParamILS: #{solved}, #{measured_runtime}, #{obj}, #{gap}, #{seed}" 
	
	#== Uncomment these lines to get access to tempfile.
	#while true
	#end
end
	
if __FILE__ == $0
	main(ARGV)
end
