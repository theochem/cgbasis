import os

__all__ = ['write_rst_table', 'write_if_changed']


def write_rst_table(f, table, nhead=1):
    """Write an RST table to file f

       **Arguments:**

       f
            A writable file-like object

       table
            A list of rows. Each row must be a list of cells containing strings.

       **Optional arguments:**

       nhead
            The number of header rows
    """

    def format_cell(cell):
        if cell is None or len(cell.strip()) == 0:
            return '\ '
        else:
            return str(cell)

    # Determine the width of each column
    widths = {}
    for row in table:
        for icell, cell in enumerate(row):
            widths[icell] = max(widths.get(icell, 2), len(format_cell(cell)))

    def format_row(cells, margin):
        return ' '.join(margin + format_cell(cell).rjust(widths[icell]) + margin
                        for icell, cell in enumerate(cells))

    # construct the column markers
    markers = format_row(['=' * widths[icell] for icell in range(len(widths))], '=')

    # top markers
    print(markers, file=f)

    # heading rows (if any)
    for irow in range(nhead):
        print(format_row(table[irow], ' '), file=f)
    if nhead > 0:
        print(markers, file=f)

    # table body
    for irow in range(nhead, len(table)):
        print(format_row(table[irow], ' '), file=f)

    # bottom markers
    print(markers, file=f)


def write_if_changed(fn, s_new):
    '''Write the string s_new to the file only if the content changes

       This is just a helper function that only really writes files when needed.
       The purpose is to avoid that Sphinx rebuilds pages that have not really
       changed.

       **Arguments:**

       fn
            The filename

       s_new
            The contents to be written to the file.
    '''
    if os.path.isfile(fn):
        # read the entire file
        with open(fn) as f:
            s_old = f.read()
        if s_new == s_old:
            print(f'File {fn} needs no update. Skipping.')
            return

    # write the new file to dis
    print(f'Writing new or updated {fn}')
    with open(fn, 'w') as f:
        f.write(s_new)
