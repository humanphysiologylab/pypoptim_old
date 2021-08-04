import time


class Timer:
    def __init__(self):
        self.times = dict()
        self.total = None

    def start(self, name) -> None:
        self.times[name] = -time.time()

    def end(self, name) -> None:
        self.times[name] += time.time()

    def report(self, sort=False) -> str:
        total = sum(self.times.values())
        self.times["total"] = total
        times = dict(
            sorted(self.times.items(), key=lambda x: x[1], reverse=True)
            if sort
            else self.times
        )
        s = "\n".join(
            f"{name}: {times[name]:.6f} {100 * times[name] / total:.2f}%"
            for name in times
        )
        del self.times["total"]
        return s

    def clear(self) -> None:
        del self.times
        self.times = dict()
