vector<string> getNextLineAndSplitIntoTokens(istream& str)
{
    vector<string>   result;
    string                line;
    getline(str,line);

    stringstream          lineStream(line);
    string                cell;

    while(getline(lineStream,cell, ' '))
    {
        result.push_back(cell);
    }
    if (!lineStream && cell.empty())
    {

        result.push_back("");
    }


    return result;
}
