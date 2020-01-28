# Function extract data from data file
# file_name must be a string which contains the extension
function extract_data(file_name)

    # Variables
    n = 0
    m = 0
    r = 0
    g = 0
    Q = 0

    # Open the text file
    f = open(file_name)

    # Read lines
    lines = readlines(f)
    for l in lines
       println(l[0])
    end

    # Close the text file
    close(f)

    # Output
    n, m, r, g, Q

end


#extract_data("E_data.txt")
