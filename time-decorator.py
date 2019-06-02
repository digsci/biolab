import time

# Time decorator from:
# https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d

def timeit(method):
  def timed(*args, **kw):
    ts = time.time()
    result = method(*args, **kw)
    te = time.time()
    print '%r  %2.2f ms' % (method.__name__, (te - ts) * 1000)
    return result
  return timed

@timeit
def f(n):
  for x in range(n):
    l = x
  return l

print(f(1000000))
